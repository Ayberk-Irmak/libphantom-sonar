// SPDX-License-Identifier: Apache-2.0
#include "phantom/ctrv.hpp"

#include <cmath>

namespace phantom {
namespace {

constexpr Real kZero = static_cast<Real>(0);
constexpr Real kOne  = static_cast<Real>(1);
constexpr Real kTwo  = static_cast<Real>(2);
constexpr Real kHugeNis = static_cast<Real>(1e30);
constexpr Real kMinRange = static_cast<Real>(1e-3);
constexpr std::size_t kN = kCtrvN;          // 5
constexpr std::size_t kMaxM = 3;

// The four quantities the coordinated turn needs, plus the two derivatives with
// respect to omega that the Jacobian needs. Computing them together is what
// keeps the series threshold in ONE place: having the transition and its
// Jacobian switch to series at different points would produce a covariance
// inconsistent with the state it describes.
struct TurnTerms {
    Real s_w;        // sin(wT)/w
    Real omc_w;      // (1 - cos(wT))/w
    Real c;          // cos(wT)
    Real s;          // sin(wT)
    Real d_s_w;      // d/dw [ sin(wT)/w ]
    Real d_omc_w;    // d/dw [ (1 - cos(wT))/w ]
};

TurnTerms turn_terms(Real w, Real dt) noexcept {
    TurnTerms t{};
    const Real wt = w * dt;
    if (std::fabs(static_cast<double>(wt)) < 1e-4) {
        const Real wt2 = wt * wt;
        const Real dt2 = dt * dt;
        const Real dt3 = dt2 * dt;
        t.s_w     = dt * (kOne - wt2 / static_cast<Real>(6));
        t.omc_w   = w * dt2 / kTwo;
        t.c       = kOne - wt2 / kTwo;
        t.s       = wt * (kOne - wt2 / static_cast<Real>(6));
        // d/dw [T - w^2 T^3/6] = -w T^3/3, and d/dw [w T^2/2] = T^2/2.
        t.d_s_w   = -w * dt3 / static_cast<Real>(3);
        t.d_omc_w = dt2 / kTwo;
    } else {
        const Real w2 = w * w;
        t.c     = std::cos(wt);
        t.s     = std::sin(wt);
        t.s_w   = t.s / w;
        t.omc_w = (kOne - t.c) / w;
        t.d_s_w   = (dt * t.c * w - t.s) / w2;
        t.d_omc_w = (dt * t.s * w - (kOne - t.c)) / w2;
    }
    return t;
}

void symmetrise(std::array<Real, kN * kN>& p) noexcept {
    for (std::size_t i = 0; i < kN; ++i) {
        for (std::size_t j = i + 1; j < kN; ++j) {
            const Real m = (p[i * kN + j] + p[j * kN + i]) / kTwo;
            p[i * kN + j] = m;
            p[j * kN + i] = m;
        }
    }
}

// Measurement Jacobian, 2 or 3 rows by 5 columns. The fifth column is zero:
// nothing measured here depends on the turn rate directly.
std::size_t jacobian(const CtrvState& s, Real r, bool with_rate,
                     std::array<Real, kMaxM * kN>& h) noexcept {
    static_assert(kMaxM * kN >= 15, "the Jacobian writes 3 rows of 5");
    for (std::size_t i = 0; i < kMaxM * kN; ++i) h[i] = kZero;
    const Real r2 = r * r;
    h[0] = s.x / r;    h[1] = s.y / r;
    h[kN + 0] = s.y / r2;  h[kN + 1] = -s.x / r2;
    if (!with_rate) return 2;
    const Real d = s.x * s.vx + s.y * s.vy;
    const Real r3 = r2 * r;
    h[2 * kN + 0] = -s.vx / r + d * s.x / r3;
    h[2 * kN + 1] = -s.vy / r + d * s.y / r3;
    h[2 * kN + 2] = -s.x / r;
    h[2 * kN + 3] = -s.y / r;
    return 3;
}

bool cholesky(std::array<Real, kMaxM * kMaxM>& a, std::size_t m) noexcept {
    for (std::size_t i = 0; i < m; ++i) {
        for (std::size_t j = 0; j <= i; ++j) {
            Real sum = a[i * m + j];
            for (std::size_t k = 0; k < j; ++k) sum -= a[i * m + k] * a[j * m + k];
            if (i == j) {
                if (!(sum > kZero)) return false;
                a[i * m + j] = std::sqrt(sum);
            } else {
                a[i * m + j] = sum / a[j * m + j];
            }
        }
        for (std::size_t j = i + 1; j < m; ++j) a[i * m + j] = kZero;
    }
    return true;
}

void cholesky_solve(const std::array<Real, kMaxM * kMaxM>& l, std::size_t m,
                    std::array<Real, kMaxM>& x) noexcept {
    for (std::size_t i = 0; i < m; ++i) {
        Real sum = x[i];
        for (std::size_t k = 0; k < i; ++k) sum -= l[i * m + k] * x[k];
        x[i] = sum / l[i * m + i];
    }
    for (std::size_t ii = m; ii-- > 0;) {
        Real sum = x[ii];
        for (std::size_t k = ii + 1; k < m; ++k) sum -= l[k * m + ii] * x[k];
        x[ii] = sum / l[ii * m + ii];
    }
}

Real range_of(const CtrvState& s) noexcept { return std::sqrt(s.x * s.x + s.y * s.y); }
Real bearing_of(const CtrvState& s) noexcept { return std::atan2(s.x, s.y); }

Real range_rate_of(const CtrvState& s) noexcept {
    const Real r = range_of(s);
    if (!(r > kMinRange)) return kZero;
    return -(s.x * s.vx + s.y * s.vy) / r;
}

std::size_t innovation(const CtrvState& s, const Measurement& z, Real r,
                       std::array<Real, kMaxM>& y) noexcept {
    y[0] = z.range_m - r;
    y[1] = wrap_pi(z.bearing_rad - bearing_of(s));
    if (!z.has_range_rate) return 2;
    y[2] = z.range_rate_mps - range_rate_of(s);
    return 3;
}

// S = H P H^T + R, and PH^T alongside it since the gain needs both.
void innovation_covariance(const CtrvTrack& t, const std::array<Real, kMaxM * kN>& h,
                           std::size_t m, const TrackerConfig& cfg,
                           std::array<Real, kN * kMaxM>& pht,
                           std::array<Real, kMaxM * kMaxM>& sm) noexcept {
    for (std::size_t i = 0; i < kN; ++i) {
        for (std::size_t j = 0; j < m; ++j) {
            Real acc = kZero;
            for (std::size_t k = 0; k < kN; ++k) acc += t.covariance[i * kN + k] * h[j * kN + k];
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
}

}  // namespace

void ctrv_initiate(CtrvTrack& track, const Measurement& z,
                   const TrackerConfig& cfg, const CtrvConfig& ct,
                   std::uint32_t id) noexcept {
    const Real sb = std::sin(z.bearing_rad);
    const Real cb = std::cos(z.bearing_rad);
    track.state.x = z.range_m * sb;
    track.state.y = z.range_m * cb;
    track.state.vx = kZero;
    track.state.vy = kZero;
    track.state.turn_rate_rps = kZero;

    for (std::size_t i = 0; i < kN * kN; ++i) track.covariance[i] = kZero;
    // Position covariance from the polar measurement, rotated into Cartesian:
    // sigma_r along the bearing, r*sigma_theta across it.
    const Real sr = cfg.range_sigma_m;
    const Real st = z.range_m * cfg.bearing_sigma_rad;
    const Real sr2 = sr * sr, st2 = st * st;
    track.covariance[0]          = sr2 * sb * sb + st2 * cb * cb;
    track.covariance[1]          = (sr2 - st2) * sb * cb;
    track.covariance[kN]         = track.covariance[1];
    track.covariance[kN + 1]     = sr2 * cb * cb + st2 * sb * sb;
    const Real sv = cfg.init_velocity_sigma_mps;
    track.covariance[2 * kN + 2] = sv * sv;
    track.covariance[3 * kN + 3] = sv * sv;
    const Real sw = ct.init_turn_rate_sigma_rps;
    track.covariance[4 * kN + 4] = sw * sw;

    track.last_update_s = z.time_s;
    track.id = id;
    track.hits = 1;
    track.misses = 0;
    track.age = 0;
    track.status = TrackStatus::Tentative;
}

void ctrv_predict(CtrvTrack& track, Real dt,
                  const TrackerConfig& cfg, const CtrvConfig& ct) noexcept {
    if (!track.live()) return;
    (void)cfg;

    CtrvState& s = track.state;
    const Real w = s.turn_rate_rps;
    const TurnTerms t = turn_terms(w, dt);

    // --- State --------------------------------------------------------------
    const Real vx = s.vx, vy = s.vy;
    s.x += t.s_w * vx - t.omc_w * vy;
    s.y += t.omc_w * vx + t.s_w * vy;
    s.vx = t.c * vx - t.s * vy;
    s.vy = t.s * vx + t.c * vy;
    // turn rate is a random walk: unchanged in the mean.

    // --- Jacobian -----------------------------------------------------------
    Real f[kN * kN];
    for (std::size_t i = 0; i < kN * kN; ++i) f[i] = kZero;
    f[0] = kOne;  f[2] = t.s_w;   f[3] = -t.omc_w;
    f[4] = t.d_s_w * vx - t.d_omc_w * vy;
    f[kN + 1] = kOne;  f[kN + 2] = t.omc_w;  f[kN + 3] = t.s_w;
    f[kN + 4] = t.d_omc_w * vx + t.d_s_w * vy;
    f[2 * kN + 2] = t.c;   f[2 * kN + 3] = -t.s;
    f[2 * kN + 4] = -dt * t.s * vx - dt * t.c * vy;
    f[3 * kN + 2] = t.s;   f[3 * kN + 3] = t.c;
    f[3 * kN + 4] = dt * t.c * vx - dt * t.s * vy;
    f[4 * kN + 4] = kOne;

    // P = F P F^T
    Real fp[kN * kN];
    for (std::size_t i = 0; i < kN; ++i) {
        for (std::size_t j = 0; j < kN; ++j) {
            Real acc = kZero;
            for (std::size_t k = 0; k < kN; ++k) acc += f[i * kN + k] * track.covariance[k * kN + j];
            fp[i * kN + j] = acc;
        }
    }
    for (std::size_t i = 0; i < kN; ++i) {
        for (std::size_t j = 0; j < kN; ++j) {
            Real acc = kZero;
            for (std::size_t k = 0; k < kN; ++k) acc += fp[i * kN + k] * f[j * kN + k];
            track.covariance[i * kN + j] = acc;
        }
    }

    // --- Q ------------------------------------------------------------------
    const Real q = ct.process_accel_mps2 * ct.process_accel_mps2;
    const Real dt2 = dt * dt, dt3 = dt2 * dt, dt4 = dt2 * dt2;
    const Real q_pp = q * dt4 / static_cast<Real>(4);
    const Real q_pv = q * dt3 / kTwo;
    const Real q_vv = q * dt2;
    track.covariance[0]           += q_pp;   track.covariance[2]           += q_pv;
    track.covariance[kN + 1]      += q_pp;   track.covariance[kN + 3]      += q_pv;
    track.covariance[2 * kN + 0]  += q_pv;   track.covariance[2 * kN + 2]  += q_vv;
    track.covariance[3 * kN + 1]  += q_pv;   track.covariance[3 * kN + 3]  += q_vv;
    // Random walk on the turn rate. This term is why omega can be learned at
    // all: without it the initial variance would shrink monotonically and the
    // filter would stop believing a manoeuvre could ever start.
    track.covariance[4 * kN + 4] +=
        ct.turn_rate_noise_rps * ct.turn_rate_noise_rps * dt;

    symmetrise(track.covariance);
    ++track.age;
}

bool ctrv_update(CtrvTrack& track, const Measurement& z,
                 const TrackerConfig& cfg, const CtrvConfig& ct) noexcept {
    if (!track.live()) return false;
    const Real r = range_of(track.state);
    if (!(r > kMinRange)) return false;

    std::array<Real, kMaxM * kN> h{};
    const std::size_t m = jacobian(track.state, r, z.has_range_rate, h);

    std::array<Real, kN * kMaxM> pht{};
    std::array<Real, kMaxM * kMaxM> sm{};
    innovation_covariance(track, h, m, cfg, pht, sm);

    std::array<Real, kMaxM * kMaxM> l = sm;
    if (!cholesky(l, m)) return false;

    Real k_gain[kN * kMaxM];
    for (std::size_t i = 0; i < kN; ++i) {
        std::array<Real, kMaxM> row{};
        for (std::size_t j = 0; j < m; ++j) row[j] = pht[i * m + j];
        cholesky_solve(l, m, row);
        for (std::size_t j = 0; j < m; ++j) k_gain[i * m + j] = row[j];
    }

    std::array<Real, kMaxM> y{};
    innovation(track.state, z, r, y);

    Real dx[kN] = {kZero, kZero, kZero, kZero, kZero};
    for (std::size_t i = 0; i < kN; ++i) {
        for (std::size_t j = 0; j < m; ++j) dx[i] += k_gain[i * m + j] * y[j];
    }
    track.state.x  += dx[0];
    track.state.y  += dx[1];
    track.state.vx += dx[2];
    track.state.vy += dx[3];
    track.state.turn_rate_rps += dx[4];

    // Clamp. The transition divides by omega, and a filter that wanders to a
    // large one predicts a tight circle, gates everything out and never
    // recovers. This is a stability guard, not a modelling choice.
    if (track.state.turn_rate_rps > ct.max_turn_rate_rps) {
        track.state.turn_rate_rps = ct.max_turn_rate_rps;
    } else if (track.state.turn_rate_rps < -ct.max_turn_rate_rps) {
        track.state.turn_rate_rps = -ct.max_turn_rate_rps;
    }

    // P = (I - K H) P
    Real kh[kN * kN];
    for (std::size_t i = 0; i < kN; ++i) {
        for (std::size_t j = 0; j < kN; ++j) {
            Real acc = kZero;
            for (std::size_t q = 0; q < m; ++q) acc += k_gain[i * m + q] * h[q * kN + j];
            kh[i * kN + j] = acc;
        }
    }
    Real np[kN * kN];
    for (std::size_t i = 0; i < kN; ++i) {
        for (std::size_t j = 0; j < kN; ++j) {
            Real acc = track.covariance[i * kN + j];
            for (std::size_t k = 0; k < kN; ++k) acc -= kh[i * kN + k] * track.covariance[k * kN + j];
            np[i * kN + j] = acc;
        }
    }
    for (std::size_t i = 0; i < kN * kN; ++i) track.covariance[i] = np[i];
    symmetrise(track.covariance);

    track.last_update_s = z.time_s;
    ++track.hits;
    track.misses = 0;
    return true;
}

Real ctrv_nis(const CtrvTrack& track, const Measurement& z,
              const TrackerConfig& cfg) noexcept {
    if (!track.live()) return kHugeNis;
    const Real r = range_of(track.state);
    if (!(r > kMinRange)) return kHugeNis;

    std::array<Real, kMaxM * kN> h{};
    const std::size_t m = jacobian(track.state, r, z.has_range_rate, h);
    std::array<Real, kN * kMaxM> pht{};
    std::array<Real, kMaxM * kMaxM> sm{};
    innovation_covariance(track, h, m, cfg, pht, sm);

    std::array<Real, kMaxM * kMaxM> l = sm;
    if (!cholesky(l, m)) return kHugeNis;

    std::array<Real, kMaxM> y{};
    innovation(track.state, z, r, y);
    std::array<Real, kMaxM> x = y;
    cholesky_solve(l, m, x);

    Real d = kZero;
    for (std::size_t i = 0; i < m; ++i) d += y[i] * x[i];
    return (d >= kZero) ? d : kHugeNis;
}

Real ctrv_turn_rate_sigma_rps(const CtrvTrack& track) noexcept {
    const Real v = track.covariance[4 * kCtrvN + 4];
    return (v > kZero) ? std::sqrt(v) : kZero;
}

}  // namespace phantom
