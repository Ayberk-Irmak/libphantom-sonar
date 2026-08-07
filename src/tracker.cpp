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

inline Real wrap_pi(Real a) noexcept {
    while (a > kPi) a -= kTwo * kPi;
    while (a < -kPi) a += kTwo * kPi;
    return a;
}

// Measurement Jacobian H (2x4), row-major, for z = [range, bearing].
void jacobian(const TargetState& s, Real r, std::array<Real, 8>& h) noexcept {
    const Real r2 = r * r;
    h[0] = s.x / r;  h[1] = s.y / r;  h[2] = kZero; h[3] = kZero;   // d range
    h[4] = s.y / r2; h[5] = -s.x / r2; h[6] = kZero; h[7] = kZero;  // d bearing
}

// S = H P H^T + R, a 2x2. Returns false if it is not invertible.
bool innovation_covariance(const Track& t, const std::array<Real, 8>& h,
                           const TrackerConfig& cfg,
                           Real& s00, Real& s01, Real& s11, Real& det) noexcept {
    // PH^T is 4x2.
    Real pht[8];
    for (std::size_t i = 0; i < kN; ++i) {
        for (std::size_t j = 0; j < 2; ++j) {
            Real acc = kZero;
            for (std::size_t k = 0; k < kN; ++k) {
                acc += t.covariance[i * kN + k] * h[j * kN + k];
            }
            pht[i * 2 + j] = acc;
        }
    }
    s00 = kZero; s01 = kZero; s11 = kZero;
    for (std::size_t k = 0; k < kN; ++k) {
        s00 += h[k] * pht[k * 2];
        s01 += h[k] * pht[k * 2 + 1];
        s11 += h[kN + k] * pht[k * 2 + 1];
    }
    s00 += cfg.range_sigma_m * cfg.range_sigma_m;
    s11 += cfg.bearing_sigma_rad * cfg.bearing_sigma_rad;

    det = s00 * s11 - s01 * s01;
    return det > kZero;
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

    std::array<Real, 8> h{};
    jacobian(track.state, r, h);

    Real s00 = 0, s01 = 0, s11 = 0, det = 0;
    if (!innovation_covariance(track, h, cfg, s00, s01, s11, det)) return kHugeNis;

    const Real y0 = z.range_m - r;
    const Real y1 = wrap_pi(z.bearing_rad - track.state.bearing_rad());

    // y^T S^-1 y with S^-1 written out for the 2x2.
    const Real inv00 = s11 / det;
    const Real inv01 = -s01 / det;
    const Real inv11 = s00 / det;
    return y0 * y0 * inv00 + kTwo * y0 * y1 * inv01 + y1 * y1 * inv11;
}

bool track_update(Track& track, const Measurement& z, const TrackerConfig& cfg) noexcept {
    if (!track.live()) return false;
    const Real r = track.state.range_m();
    if (!(r > kMinRange)) return false;

    std::array<Real, 8> h{};
    jacobian(track.state, r, h);

    Real s00 = 0, s01 = 0, s11 = 0, det = 0;
    if (!innovation_covariance(track, h, cfg, s00, s01, s11, det)) return false;

    const Real inv00 = s11 / det;
    const Real inv01 = -s01 / det;
    const Real inv11 = s00 / det;

    // K = P H^T S^-1, 4x2.
    Real pht[8];
    for (std::size_t i = 0; i < kN; ++i) {
        for (std::size_t j = 0; j < 2; ++j) {
            Real acc = kZero;
            for (std::size_t k = 0; k < kN; ++k) {
                acc += track.covariance[i * kN + k] * h[j * kN + k];
            }
            pht[i * 2 + j] = acc;
        }
    }
    Real k_gain[8];
    for (std::size_t i = 0; i < kN; ++i) {
        k_gain[i * 2]     = pht[i * 2] * inv00 + pht[i * 2 + 1] * inv01;
        k_gain[i * 2 + 1] = pht[i * 2] * inv01 + pht[i * 2 + 1] * inv11;
    }

    const Real y0 = z.range_m - r;
    const Real y1 = wrap_pi(z.bearing_rad - track.state.bearing_rad());

    track.state.x  += k_gain[0] * y0 + k_gain[1] * y1;
    track.state.y  += k_gain[2] * y0 + k_gain[3] * y1;
    track.state.vx += k_gain[4] * y0 + k_gain[5] * y1;
    track.state.vy += k_gain[6] * y0 + k_gain[7] * y1;

    // P = (I - K H) P
    Real kh[16];
    for (std::size_t i = 0; i < kN; ++i) {
        for (std::size_t j = 0; j < kN; ++j) {
            kh[i * kN + j] = k_gain[i * 2] * h[j] + k_gain[i * 2 + 1] * h[kN + j];
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

    // --- 2. Greedy nearest-neighbour association ----------------------------
    // At most 64 measurements are considered per step; beyond that the extras
    // are treated as unassociated and may start tracks, which is the same
    // outcome a global assignment would give for a cluttered scan.
    constexpr std::size_t kMaxMeas = 64;
    bool used[kMaxMeas] = {};
    const std::size_t n_meas = (measurements.size() < kMaxMeas) ? measurements.size() : kMaxMeas;

    for (Track& t : tracks) {
        if (!t.live()) continue;
        std::size_t best = n_meas;
        Real best_nis = cfg.gate_chi2;
        for (std::size_t m = 0; m < n_meas; ++m) {
            if (used[m]) continue;
            const Real d = track_nis(t, measurements[m], cfg);
            if (d < best_nis) { best_nis = d; best = m; }
        }
        if (best < n_meas && track_update(t, measurements[best], cfg)) {
            used[best] = true;
        } else {
            ++t.misses;
        }
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

std::size_t count_tracks(std::span<const Track> tracks, TrackStatus status) noexcept {
    std::size_t n = 0;
    for (const Track& t : tracks) {
        if (t.status == status) ++n;
    }
    return n;
}

}  // namespace phantom
