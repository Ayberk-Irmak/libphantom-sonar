// SPDX-License-Identifier: Apache-2.0
#include "phantom/tracker.hpp"

#include <cmath>

namespace phantom {
namespace {

constexpr Real kZero = static_cast<Real>(0);
constexpr Real kOne  = static_cast<Real>(1);
constexpr Real kTwo  = static_cast<Real>(2);

constexpr std::size_t kN = 4;   // state dimension
constexpr Real kHugeNis = static_cast<Real>(1e12);

// A track closer to the origin than this has no defined bearing, and the
// measurement Jacobian divides by r^2.
constexpr Real kMinRange = static_cast<Real>(1e-3);

constexpr std::size_t kMaxM = 3;   // range, bearing, and optionally range rate

// The measurement Jacobian is kMaxM rows of kN. Named rather than written as
// the product at each use: cppcheck 2.13 (what CI ships) resolves
// std::array<Real, kJacobianSize> as having THREE elements -- it takes the first
// factor and stops -- and then reports every write past index 2 as out of
// bounds. The code is correct at 12 elements and gcc, clang and ASan all agree,
// but a name the analyser can fold is cheaper than arguing with it.
constexpr std::size_t kJacobianSize = 12;
constexpr std::size_t kInnovSize = 9;
constexpr std::size_t kPhtSize = 12;
static_assert(kPhtSize == kN * kMaxM, "P H^T is n x m");
static_assert(kInnovSize == kMaxM * kMaxM, "the innovation covariance is m x m");
static_assert(kJacobianSize == kMaxM * kN, "3 measurement rows by 4 states");

// Measurement Jacobian H (m x 4), row-major, for z = [range, bearing, rdot].
//
// The third row is where the fusion lives. With rdot = -(x vx + y vy)/r:
//   d/dx  = -vx/r + d x / r^3      d/dvx = -x/r
//   d/dy  = -vy/r + d y / r^3      d/dvy = -y/r
// It is the only row that touches the velocity states, which is exactly why a
// direct range-rate measurement observes velocity that position history can
// only infer.
std::size_t jacobian(const TargetState& s, Real r, bool with_rate,
                     std::array<Real, kJacobianSize>& h) noexcept {
    const Real r2 = r * r;
    h[0] = s.x / r;   h[1] = s.y / r;   h[2] = kZero; h[3] = kZero;
    h[4] = s.y / r2;  h[5] = -s.x / r2; h[6] = kZero; h[7] = kZero;
    if (!with_rate) return 2;
    const Real d = s.x * s.vx + s.y * s.vy;
    const Real r3 = r2 * r;
    h[8]  = -s.vx / r + d * s.x / r3;
    h[9]  = -s.vy / r + d * s.y / r3;
    h[10] = -s.x / r;
    h[11] = -s.y / r;
    return 3;
}

// Cholesky factor of an m x m symmetric positive definite matrix, in place over
// the lower triangle. Small and real, so no need for the complex version in
// beamformer.cpp.
bool small_cholesky(std::array<Real, kInnovSize>& a, std::size_t m) noexcept {
    for (std::size_t j = 0; j < m; ++j) {
        Real d = a[j * m + j];
        for (std::size_t k = 0; k < j; ++k) d -= a[j * m + k] * a[j * m + k];
        if (!(d > kZero)) return false;
        const Real ljj = std::sqrt(d);
        a[j * m + j] = ljj;
        for (std::size_t i = j + 1; i < m; ++i) {
            Real acc = a[i * m + j];
            for (std::size_t k = 0; k < j; ++k) acc -= a[i * m + k] * a[j * m + k];
            a[i * m + j] = acc / ljj;
        }
    }
    return true;
}

void small_cholesky_solve(const std::array<Real, kInnovSize>& l, std::size_t m,
                          std::array<Real, kMaxM>& x) noexcept {
    for (std::size_t i = 0; i < m; ++i) {
        Real acc = x[i];
        for (std::size_t k = 0; k < i; ++k) acc -= l[i * m + k] * x[k];
        x[i] = acc / l[i * m + i];
    }
    for (std::size_t ii = m; ii-- > 0;) {
        Real acc = x[ii];
        for (std::size_t k = ii + 1; k < m; ++k) acc -= l[k * m + ii] * x[k];
        x[ii] = acc / l[ii * m + ii];
    }
}

// S = H P H^T + R and the PH^T it is built from, both needed by the update.
bool innovation_covariance(const Track& t, const std::array<Real, kJacobianSize>& h,
                           std::size_t m, const TrackerConfig& cfg,
                           std::array<Real, kPhtSize>& pht,
                           std::array<Real, kInnovSize>& sm) noexcept {
    for (std::size_t i = 0; i < kN; ++i) {
        for (std::size_t j = 0; j < m; ++j) {
            Real acc = kZero;
            for (std::size_t k = 0; k < kN; ++k) {
                acc += t.covariance[i * kN + k] * h[j * kN + k];
            }
            pht[i * m + j] = acc;
        }
    }
    for (std::size_t i = 0; i < m; ++i) {
        for (std::size_t j = 0; j < m; ++j) {
            Real acc = kZero;
            for (std::size_t k = 0; k < kN; ++k) acc += h[i * kN + k] * pht[k * m + j];
            sm[i * m + j] = acc;
        }
    }
    // Add R to the diagonal. Written as a loop over the diagonal rather than
    // three indexed writes: sm[0], sm[m+1], sm[8] is correct for m = 2 and 3 but
    // only provably so if you know m cannot exceed 3, which a reader -- and a
    // static analyser -- has to take on trust. This form is obviously in bounds.
    const Real r_diag[kMaxM] = {cfg.range_sigma_m * cfg.range_sigma_m,
                                cfg.bearing_sigma_rad * cfg.bearing_sigma_rad,
                                cfg.range_rate_sigma_mps * cfg.range_rate_sigma_mps};
    for (std::size_t i = 0; i < m && i < kMaxM; ++i) sm[i * m + i] += r_diag[i];
    return true;
}

// Innovation vector, with the bearing wrapped.
std::size_t innovation(const Track& t, const Measurement& z, Real r,
                       std::array<Real, kMaxM>& y) noexcept {
    y[0] = z.range_m - r;
    y[1] = wrap_pi(z.bearing_rad - t.state.bearing_rad());
    if (!z.has_range_rate) return 2;
    y[2] = z.range_rate_mps - t.state.range_rate_mps();
    return 3;
}

}  // namespace

Real TargetState::range_m() const noexcept { return std::sqrt(x * x + y * y); }

Real TargetState::bearing_rad() const noexcept { return std::atan2(x, y); }

Real TargetState::range_rate_mps() const noexcept {
    const Real r = range_m();
    if (!(r > kMinRange)) return kZero;
    // dr/dt = (x vx + y vy) / r; closing is a shrinking range, hence the sign.
    return -(x * vx + y * vy) / r;
}

Real chi2_gate_2dof(Real probability) noexcept {
    if (!(probability > kZero) || !(probability < kOne)) return kZero;
    // The chi-square CDF with 2 dof is 1 - exp(-x/2), so it inverts in closed
    // form: no table, no approximation.
    return -kTwo * std::log(kOne - probability);
}

Real chi2_gate(Real probability, std::size_t dof) noexcept {
    if (!(probability > kZero) || !(probability < kOne) || dof == 0) return kZero;
    if (dof == 2) return chi2_gate_2dof(probability);
    // No closed form for odd dof, so bisect on the CDF. The 3-dof CDF is
    // erf(sqrt(x/2)) - sqrt(2x/pi) exp(-x/2); for 1 dof it is erf(sqrt(x/2)).
    auto cdf = [dof](Real x) noexcept -> Real {
        const Real a = std::erf(std::sqrt(x / kTwo));
        if (dof == 1) return a;
        if (dof == 3) {
            return a - std::sqrt(kTwo * x / kPi) * std::exp(-x / kTwo);
        }
        return a;   // higher dof are not needed here
    };
    Real lo = kZero;
    Real hi = static_cast<Real>(100);
    for (int i = 0; i < 100; ++i) {
        const Real mid = (lo + hi) / kTwo;
        if (cdf(mid) < probability) lo = mid; else hi = mid;
    }
    return lo;
}

void track_predict(Track& track, Real dt, const TrackerConfig& cfg) noexcept {
    if (!track.live()) return;

    track.state.x += track.state.vx * dt;
    track.state.y += track.state.vy * dt;

    // P = F P F^T + Q with F the constant-velocity transition. Written out
    // rather than looped: F is mostly identity, and the explicit form makes the
    // structure checkable against the textbook.
    std::array<Real, 16>& p = track.covariance;
    // F P
    Real fp[16];
    for (std::size_t j = 0; j < kN; ++j) {
        fp[0 * kN + j] = p[0 * kN + j] + dt * p[2 * kN + j];
        fp[1 * kN + j] = p[1 * kN + j] + dt * p[3 * kN + j];
        fp[2 * kN + j] = p[2 * kN + j];
        fp[3 * kN + j] = p[3 * kN + j];
    }
    // (F P) F^T
    for (std::size_t i = 0; i < kN; ++i) {
        p[i * kN + 0] = fp[i * kN + 0] + dt * fp[i * kN + 2];
        p[i * kN + 1] = fp[i * kN + 1] + dt * fp[i * kN + 3];
        p[i * kN + 2] = fp[i * kN + 2];
        p[i * kN + 3] = fp[i * kN + 3];
    }

    // Discrete white-noise acceleration.
    const Real q = cfg.process_accel_mps2 * cfg.process_accel_mps2;
    const Real dt2 = dt * dt;
    const Real dt3 = dt2 * dt;
    const Real dt4 = dt2 * dt2;
    const Real q_pp = q * dt4 / static_cast<Real>(4);
    const Real q_pv = q * dt3 / kTwo;
    const Real q_vv = q * dt2;

    p[0 * kN + 0] += q_pp;  p[0 * kN + 2] += q_pv;
    p[1 * kN + 1] += q_pp;  p[1 * kN + 3] += q_pv;
    p[2 * kN + 0] += q_pv;  p[2 * kN + 2] += q_vv;
    p[3 * kN + 1] += q_pv;  p[3 * kN + 3] += q_vv;

    ++track.age;
}

Real track_nis(const Track& track, const Measurement& z,
               const TrackerConfig& cfg) noexcept {
    if (!track.live()) return kHugeNis;
    const Real r = track.state.range_m();
    if (!(r > kMinRange)) return kHugeNis;

    std::array<Real, kJacobianSize> h{};
    const std::size_t m = jacobian(track.state, r, z.has_range_rate, h);

    std::array<Real, kPhtSize> pht{};
    std::array<Real, kInnovSize> sm{};
    innovation_covariance(track, h, m, cfg, pht, sm);

    std::array<Real, kMaxM> y{};
    innovation(track, z, r, y);

    std::array<Real, kInnovSize> l = sm;
    if (!small_cholesky(l, m)) return kHugeNis;
    std::array<Real, kMaxM> x = y;
    small_cholesky_solve(l, m, x);

    Real d = kZero;
    for (std::size_t i = 0; i < m; ++i) d += y[i] * x[i];
    return (d >= kZero) ? d : kHugeNis;
}

Real measurement_likelihood(const Track& track, const Measurement& z,
                            const TrackerConfig& cfg) noexcept {
    if (!track.live()) return kZero;
    const Real r = track.state.range_m();
    if (!(r > kMinRange)) return kZero;

    std::array<Real, kJacobianSize> h{};
    const std::size_t m = jacobian(track.state, r, z.has_range_rate, h);

    std::array<Real, kPhtSize> pht{};
    std::array<Real, kInnovSize> sm{};
    innovation_covariance(track, h, m, cfg, pht, sm);

    std::array<Real, kInnovSize> l = sm;
    if (!small_cholesky(l, m)) return kZero;

    std::array<Real, kMaxM> y{};
    innovation(track, z, r, y);
    std::array<Real, kMaxM> x = y;
    small_cholesky_solve(l, m, x);

    Real d = kZero;
    for (std::size_t i = 0; i < m; ++i) d += y[i] * x[i];
    if (!(d >= kZero)) return kZero;

    // |S| = (prod L_ii)^2 from the Cholesky factor, so the determinant costs
    // nothing beyond the solve that was needed anyway.
    Real log_det = kZero;
    for (std::size_t i = 0; i < m; ++i) {
        const Real lii = l[i * m + i];
        if (!(lii > kZero)) return kZero;
        log_det += kTwo * std::log(lii);
    }
    const Real log_two_pi = static_cast<Real>(1.8378770664093454836);
    const Real exponent = -(d + log_det + static_cast<Real>(m) * log_two_pi) / kTwo;
    // Guard the exponential: a model that fits terribly must return 0, not a
    // denormal that later divides into something.
    if (exponent < static_cast<Real>(-700)) return kZero;
    return std::exp(exponent);
}

bool track_update(Track& track, const Measurement& z, const TrackerConfig& cfg) noexcept {
    if (!track.live()) return false;
    const Real r = track.state.range_m();
    if (!(r > kMinRange)) return false;

    std::array<Real, kJacobianSize> h{};
    const std::size_t m = jacobian(track.state, r, z.has_range_rate, h);

    std::array<Real, kPhtSize> pht{};
    std::array<Real, kInnovSize> sm{};
    innovation_covariance(track, h, m, cfg, pht, sm);

    std::array<Real, kInnovSize> l = sm;
    if (!small_cholesky(l, m)) return false;

    // K = P H^T S^-1, obtained a column at a time by solving S k_i = (PH^T)_i.
    Real k_gain[kN * kMaxM];
    for (std::size_t i = 0; i < kN; ++i) {
        std::array<Real, kMaxM> row{};
        for (std::size_t j = 0; j < m; ++j) row[j] = pht[i * m + j];
        small_cholesky_solve(l, m, row);
        for (std::size_t j = 0; j < m; ++j) k_gain[i * m + j] = row[j];
    }

    std::array<Real, kMaxM> y{};
    innovation(track, z, r, y);

    Real dx[kN] = {kZero, kZero, kZero, kZero};
    for (std::size_t i = 0; i < kN; ++i) {
        for (std::size_t j = 0; j < m; ++j) dx[i] += k_gain[i * m + j] * y[j];
    }
    track.state.x  += dx[0];
    track.state.y  += dx[1];
    track.state.vx += dx[2];
    track.state.vy += dx[3];

    // P = (I - K H) P
    Real kh[16];
    for (std::size_t i = 0; i < kN; ++i) {
        for (std::size_t j = 0; j < kN; ++j) {
            Real acc = kZero;
            for (std::size_t q = 0; q < m; ++q) acc += k_gain[i * m + q] * h[q * kN + j];
            kh[i * kN + j] = acc;
        }
    }
    Real np[16];
    for (std::size_t i = 0; i < kN; ++i) {
        for (std::size_t j = 0; j < kN; ++j) {
            Real acc = track.covariance[i * kN + j];
            for (std::size_t k = 0; k < kN; ++k) {
                acc -= kh[i * kN + k] * track.covariance[k * kN + j];
            }
            np[i * kN + j] = acc;
        }
    }
    // Symmetrise. (I-KH)P is symmetric in exact arithmetic and drifts out of it
    // in floating point; an asymmetric covariance eventually goes indefinite
    // and the filter diverges with no warning.
    for (std::size_t i = 0; i < kN; ++i) {
        for (std::size_t j = 0; j < kN; ++j) {
            track.covariance[i * kN + j] = (np[i * kN + j] + np[j * kN + i]) / kTwo;
        }
    }

    track.last_update_s = z.time_s;
    // Exponentially smoothed amplitude, for ghost recognition.
    track.amplitude = (track.hits == 0)
                    ? z.amplitude
                    : static_cast<Real>(0.7) * track.amplitude
                    + static_cast<Real>(0.3) * z.amplitude;
    track.label = z.label;
    ++track.hits;
    track.misses = 0;
    return true;
}

void track_initiate(Track& track, const Measurement& z, const TrackerConfig& cfg,
                    std::uint32_t id) noexcept {
    track.state.x = z.range_m * std::sin(z.bearing_rad);
    track.state.y = z.range_m * std::cos(z.bearing_rad);
    track.state.vx = kZero;
    track.state.vy = kZero;

    for (Real& v : track.covariance) v = kZero;
    // Position uncertainty from the measurement, converted to Cartesian: the
    // cross-range error is r*sigma_bearing, which at long range dwarfs the
    // range error and is why bearing accuracy matters so much.
    const Real cross = z.range_m * cfg.bearing_sigma_rad;
    const Real rad = cfg.range_sigma_m;
    const Real var_pos = (cross * cross + rad * rad);   // isotropic upper bound
    track.covariance[0] = var_pos;
    track.covariance[5] = var_pos;
    track.covariance[10] = cfg.init_velocity_sigma_mps * cfg.init_velocity_sigma_mps;
    track.covariance[15] = cfg.init_velocity_sigma_mps * cfg.init_velocity_sigma_mps;

    track.last_update_s = z.time_s;
    track.amplitude = z.amplitude;
    track.label = z.label;
    track.id = id;
    track.hits = 1;
    track.misses = 0;
    track.age = 0;
    track.status = TrackStatus::Tentative;
}

std::size_t tracker_step(std::span<Track> tracks,
                         std::span<const Measurement> measurements,
                         const TrackerConfig& cfg,
                         Real time_s,
                         std::uint32_t& next_id) noexcept {
    // --- 1. Predict every live track to now ---------------------------------
    for (Track& t : tracks) {
        if (!t.live()) continue;
        const Real dt = time_s - t.last_update_s;
        if (dt > kZero) track_predict(t, dt, cfg);
    }

    // --- 2. Association, in global cost order -------------------------------
    // Per-track greedy assigns in TRACK order, so the first track takes the
    // measurement it likes best even when a later track wants it far more --
    // which is exactly how two crossing targets swap identities. Ordering the
    // whole (track, measurement) set by NIS and assigning best-first removes
    // that particular failure.
    //
    // Still not the optimal assignment: a case exists where the globally
    // cheapest pair forces two expensive ones and a Hungarian solve would do
    // better. But it is O(T*M log(T*M)) with no allocation, and it fixes the
    // crossing swap that per-track greedy produces routinely.
    constexpr std::size_t kMaxMeas = 64;
    constexpr std::size_t kMaxPairs = 512;
    bool used[kMaxMeas] = {};
    const std::size_t n_meas = (measurements.size() < kMaxMeas) ? measurements.size() : kMaxMeas;

    struct Pair { Real cost; std::uint16_t track; std::uint16_t meas; };
    Pair pairs[kMaxPairs];
    std::size_t n_pairs = 0;

    for (std::size_t ti = 0; ti < tracks.size() && n_pairs < kMaxPairs; ++ti) {
        if (!tracks[ti].live()) continue;
        for (std::size_t m = 0; m < n_meas && n_pairs < kMaxPairs; ++m) {
            const Real d = track_nis(tracks[ti], measurements[m], cfg);
            const Real gate = measurements[m].has_range_rate ? cfg.gate_chi2_3dof
                                                             : cfg.gate_chi2_2dof;
            if (d < gate) {
                pairs[n_pairs++] = Pair{d, static_cast<std::uint16_t>(ti),
                                        static_cast<std::uint16_t>(m)};
            }
        }
    }

    // Insertion sort: n_pairs is small and bounded, and this needs no scratch.
    for (std::size_t i = 1; i < n_pairs; ++i) {
        const Pair key = pairs[i];
        std::size_t j = i;
        while (j > 0 && pairs[j - 1].cost > key.cost) {
            pairs[j] = pairs[j - 1];
            --j;
        }
        pairs[j] = key;
    }

    bool track_taken[kMaxMeas] = {};
    const std::size_t n_tracks = (tracks.size() < kMaxMeas) ? tracks.size() : kMaxMeas;
    for (std::size_t i = 0; i < n_pairs; ++i) {
        const std::size_t ti = pairs[i].track;
        const std::size_t mi = pairs[i].meas;
        if (ti >= n_tracks || track_taken[ti] || used[mi]) continue;
        if (track_update(tracks[ti], measurements[mi], cfg)) {
            track_taken[ti] = true;
            used[mi] = true;
        }
    }
    for (std::size_t ti = 0; ti < tracks.size(); ++ti) {
        if (!tracks[ti].live()) continue;
        if (ti >= n_tracks || !track_taken[ti]) ++tracks[ti].misses;
    }

    // --- 3. Confirmation and deletion ---------------------------------------
    for (Track& t : tracks) {
        if (!t.live()) continue;
        if (t.misses >= cfg.delete_misses) {
            t.status = TrackStatus::Free;
            continue;
        }
        if (t.status == TrackStatus::Tentative && t.hits >= cfg.confirm_hits) {
            t.status = TrackStatus::Confirmed;
        } else if (t.misses > 0 && t.status == TrackStatus::Confirmed) {
            t.status = TrackStatus::Coasting;
        } else if (t.misses == 0 && t.status == TrackStatus::Coasting) {
            t.status = TrackStatus::Confirmed;
        }
    }

    // --- 4. Initiate from what nothing claimed ------------------------------
    for (std::size_t m = 0; m < n_meas; ++m) {
        if (used[m]) continue;
        for (Track& t : tracks) {
            if (t.live()) continue;
            track_initiate(t, measurements[m], cfg, next_id++);
            break;
        }
    }

    std::size_t live = 0;
    for (const Track& t : tracks) {
        if (t.live()) ++live;
    }
    return live;
}

std::size_t suppress_template_ghosts(std::span<Track> tracks,
                                     const GhostConfig& cfg) noexcept {
    std::size_t killed = 0;
    for (std::size_t i = 0; i < tracks.size(); ++i) {
        Track& a = tracks[i];
        if (!a.live() || a.hits < cfg.min_hits) continue;
        for (std::size_t j = i + 1; j < tracks.size(); ++j) {
            Track& b = tracks[j];
            if (!b.live() || b.hits < cfg.min_hits) continue;

            // THE safety check. A ghost is the same arrival seen through a
            // different matched filter, so the labels must differ. Two real
            // targets illuminated by one sonar return the same waveform and are
            // therefore never paired -- without this, a formation in line
            // astern would be deleted, and the test suite demonstrates that.
            if (a.label == b.label) continue;

            const Real dbear = std::fabs(wrap_pi(a.state.bearing_rad()
                                               - b.state.bearing_rad()));
            if (dbear > cfg.max_bearing_delta_rad) continue;

            const Real drate = std::fabs(a.state.range_rate_mps()
                                       - b.state.range_rate_mps());
            if (drate > cfg.max_range_rate_delta_mps) continue;

            const Real doffset = std::fabs(a.state.range_m() - b.state.range_m());
            if (!(doffset > kZero) || doffset > cfg.max_range_offset_m) continue;

            Track* weak = (a.amplitude < b.amplitude) ? &a : &b;
            const Track* strong = (weak == &a) ? &b : &a;
            if (!(weak->amplitude > kZero) || !(strong->amplitude > kZero)) continue;
            const Real ratio_db = static_cast<Real>(20)
                                * std::log10(strong->amplitude / weak->amplitude);
            if (ratio_db < cfg.min_amplitude_ratio_db) continue;

            weak->status = TrackStatus::Free;
            ++killed;
            if (weak == &a) break;   // `a` is gone; move to the next i
        }
    }
    return killed;
}

std::size_t count_tracks(std::span<const Track> tracks, TrackStatus status) noexcept {
    std::size_t n = 0;
    for (const Track& t : tracks) {
        if (t.status == status) ++n;
    }
    return n;
}

std::size_t count_established(std::span<const Track> tracks) noexcept {
    std::size_t n = 0;
    for (const Track& t : tracks) {
        if (t.status == TrackStatus::Confirmed || t.status == TrackStatus::Coasting) ++n;
    }
    return n;
}

}  // namespace phantom
