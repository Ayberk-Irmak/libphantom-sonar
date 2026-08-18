// SPDX-License-Identifier: Apache-2.0
#include "phantom/imm.hpp"

#include <cmath>

namespace phantom {
namespace {

constexpr Real kZero = static_cast<Real>(0);
constexpr Real kOne  = static_cast<Real>(1);
constexpr Real kTwo  = static_cast<Real>(2);
constexpr std::size_t kN = 4;
// Matches tracker.cpp's: a NIS no real measurement can reach, so it never gates.
constexpr Real kHugeNis = static_cast<Real>(1e30);

// Turn rate of each model, in units of imm.turn_rate_rps.
constexpr Real kModelRateScale[kImmModels] = {kZero, kOne, -kOne};

Real* model_cov(ImmTrack& t, std::size_t j) noexcept { return &t.model_cov[j * 16]; }

void state_to_array(const TargetState& s, Real (&v)[kN]) noexcept {
    v[0] = s.x; v[1] = s.y; v[2] = s.vx; v[3] = s.vy;
}

void array_to_state(const Real (&v)[kN], TargetState& s) noexcept {
    s.x = v[0]; s.y = v[1]; s.vx = v[2]; s.vy = v[3];
}

// Coordinated-turn transition for turn rate w over dt. At w = 0 this is exactly
// the constant-velocity matrix, which is why the CV model is not a special
// case in the code -- but the limit has to be taken by hand, because sin(w dt)/w
// is 0/0 there and evaluating it in floating point loses all precision long
// before w reaches zero.
void transition(Real w, Real dt, Real (&f)[16]) noexcept {
    for (std::size_t i = 0; i < 16; ++i) f[i] = kZero;
    f[0] = kOne; f[5] = kOne;

    Real s_over_w, one_minus_c_over_w, cw, sw;
    // The threshold is where the series and the direct form agree to working
    // precision. Below it the series is not an approximation to the answer, it
    // is the more accurate way to compute it.
    if (std::fabs(static_cast<double>(w * dt)) < 1e-4) {
        const Real wt = w * dt;
        const Real wt2 = wt * wt;
        s_over_w           = dt * (kOne - wt2 / static_cast<Real>(6));
        one_minus_c_over_w = dt * wt / kTwo;
        cw = kOne - wt2 / kTwo;
        sw = wt * (kOne - wt2 / static_cast<Real>(6));
    } else {
        cw = std::cos(w * dt);
        sw = std::sin(w * dt);
        s_over_w           = sw / w;
        one_minus_c_over_w = (kOne - cw) / w;
    }

    f[2]  =  s_over_w;            f[3]  = -one_minus_c_over_w;
    f[6]  =  one_minus_c_over_w;  f[7]  =  s_over_w;
    f[10] =  cw;                  f[11] = -sw;
    f[14] =  sw;                  f[15] =  cw;
}

// P <- F P F^T + Q, with Q the discrete white-noise acceleration for `sigma_a`.
void propagate(const Real (&f)[16], Real sigma_a, Real dt, Real* p) noexcept {
    Real fp[16];
    for (std::size_t i = 0; i < kN; ++i) {
        for (std::size_t j = 0; j < kN; ++j) {
            Real acc = kZero;
            for (std::size_t k = 0; k < kN; ++k) acc += f[i * kN + k] * p[k * kN + j];
            fp[i * kN + j] = acc;
        }
    }
    for (std::size_t i = 0; i < kN; ++i) {
        for (std::size_t j = 0; j < kN; ++j) {
            Real acc = kZero;
            for (std::size_t k = 0; k < kN; ++k) acc += fp[i * kN + k] * f[j * kN + k];
            p[i * kN + j] = acc;
        }
    }
    const Real q = sigma_a * sigma_a;
    const Real dt2 = dt * dt, dt3 = dt2 * dt, dt4 = dt2 * dt2;
    const Real q_pp = q * dt4 / static_cast<Real>(4);
    const Real q_pv = q * dt3 / kTwo;
    const Real q_vv = q * dt2;
    p[0]  += q_pp;  p[2]  += q_pv;
    p[5]  += q_pp;  p[7]  += q_pv;
    p[8]  += q_pv;  p[10] += q_vv;
    p[13] += q_pv;  p[15] += q_vv;
    // Symmetrise, for the same reason track_update does: an asymmetric
    // covariance eventually goes indefinite and the filter diverges silently.
    for (std::size_t i = 0; i < kN; ++i) {
        for (std::size_t j = i + 1; j < kN; ++j) {
            const Real m = (p[i * kN + j] + p[j * kN + i]) / kTwo;
            p[i * kN + j] = m;
            p[j * kN + i] = m;
        }
    }
}

Real model_sigma(const ImmConfig& imm, std::size_t j) noexcept {
    return (j == 0) ? imm.cv_accel_mps2 : imm.ct_accel_mps2;
}

// Markov transition probability from model i to model j.
Real switch_prob(const ImmConfig& imm, std::size_t i, std::size_t j) noexcept {
    const Real off = imm.switch_probability / static_cast<Real>(kImmModels - 1);
    return (i == j) ? (kOne - imm.switch_probability) : off;
}

}  // namespace

void imm_initiate(ImmTrack& track, const Measurement& z,
                  const TrackerConfig& cfg, const ImmConfig& imm,
                  std::uint32_t id) noexcept {
    Track seed;
    track_initiate(seed, z, cfg, id);
    track.state = seed.state;
    track.covariance = seed.covariance;
    for (std::size_t j = 0; j < kImmModels; ++j) {
        track.model_state[j] = seed.state;
        for (std::size_t k = 0; k < 16; ++k) model_cov(track, j)[k] = seed.covariance[k];
    }
    // A target is presumed straight until it shows otherwise. Splitting the
    // remainder evenly between the two turn directions is the only neutral
    // choice: nothing in one detection says which way it might go.
    track.model_prob[0] = static_cast<Real>(0.8);
    track.model_prob[1] = static_cast<Real>(0.1);
    track.model_prob[2] = static_cast<Real>(0.1);
    track.last_update_s = z.time_s;
    track.amplitude = z.amplitude;
    track.label = z.label;
    track.id = id;
    track.hits = 1;
    track.misses = 0;
    track.age = 0;
    track.status = TrackStatus::Tentative;
    (void)imm;
}

void imm_predict(ImmTrack& track, Real dt,
                 const TrackerConfig& cfg, const ImmConfig& imm) noexcept {
    if (!track.live()) return;
    (void)cfg;

    // --- Mixing -----------------------------------------------------------
    // c_j = sum_i pi_ij mu_i, and mu_{i|j} = pi_ij mu_i / c_j.
    Real cbar[kImmModels];
    for (std::size_t j = 0; j < kImmModels; ++j) {
        Real acc = kZero;
        for (std::size_t i = 0; i < kImmModels; ++i) {
            acc += switch_prob(imm, i, j) * track.model_prob[i];
        }
        cbar[j] = acc;
    }

    TargetState mixed_state[kImmModels];
    Real mixed_cov[kImmModels][16];
    for (std::size_t j = 0; j < kImmModels; ++j) {
        Real xm[kN] = {kZero, kZero, kZero, kZero};
        Real w[kImmModels];
        for (std::size_t i = 0; i < kImmModels; ++i) {
            w[i] = (cbar[j] > kZero)
                 ? switch_prob(imm, i, j) * track.model_prob[i] / cbar[j]
                 : (i == j ? kOne : kZero);
            Real xi[kN];
            state_to_array(track.model_state[i], xi);
            for (std::size_t k = 0; k < kN; ++k) xm[k] += w[i] * xi[k];
        }
        array_to_state(xm, mixed_state[j]);

        // P0_j = sum_i mu_{i|j} [ P_i + (x_i - x0_j)(x_i - x0_j)^T ]. The outer
        // product is the part that matters: without it the mixed covariance
        // ignores how far apart the models' estimates are, and the IMM becomes
        // overconfident exactly when the models disagree -- which is exactly
        // when it should not be.
        for (std::size_t k = 0; k < 16; ++k) mixed_cov[j][k] = kZero;
        for (std::size_t i = 0; i < kImmModels; ++i) {
            if (w[i] <= kZero) continue;
            Real xi[kN];
            state_to_array(track.model_state[i], xi);
            Real d[kN];
            for (std::size_t k = 0; k < kN; ++k) d[k] = xi[k] - xm[k];
            const Real* pi = &track.model_cov[i * 16];
            for (std::size_t a = 0; a < kN; ++a) {
                for (std::size_t b = 0; b < kN; ++b) {
                    mixed_cov[j][a * kN + b] += w[i] * (pi[a * kN + b] + d[a] * d[b]);
                }
            }
        }
    }

    // --- Per-model prediction ---------------------------------------------
    for (std::size_t j = 0; j < kImmModels; ++j) {
        Real f[16];
        transition(kModelRateScale[j] * imm.turn_rate_rps, dt, f);

        Real x[kN];
        state_to_array(mixed_state[j], x);
        Real xn[kN] = {kZero, kZero, kZero, kZero};
        for (std::size_t a = 0; a < kN; ++a) {
            for (std::size_t b = 0; b < kN; ++b) xn[a] += f[a * kN + b] * x[b];
        }
        array_to_state(xn, track.model_state[j]);

        Real* p = model_cov(track, j);
        for (std::size_t k = 0; k < 16; ++k) p[k] = mixed_cov[j][k];
        propagate(f, model_sigma(imm, j), dt, p);
    }

    // The combined prediction, so a caller can gate before updating.
    imm_update_combined_estimate(track);
    ++track.age;
}

bool imm_update(ImmTrack& track, const Measurement& z,
                const TrackerConfig& cfg, const ImmConfig& imm) noexcept {
    if (!track.live()) return false;

    // c_j again: the probability of being in model j after the switch, before
    // the measurement is seen. This is the prior the likelihood multiplies.
    Real cbar[kImmModels];
    for (std::size_t j = 0; j < kImmModels; ++j) {
        Real acc = kZero;
        for (std::size_t i = 0; i < kImmModels; ++i) {
            acc += switch_prob(imm, i, j) * track.model_prob[i];
        }
        cbar[j] = acc;
    }

    Real like[kImmModels];
    bool any = false;
    for (std::size_t j = 0; j < kImmModels; ++j) {
        Track tmp;
        tmp.state = track.model_state[j];
        for (std::size_t k = 0; k < 16; ++k) tmp.covariance[k] = model_cov(track, j)[k];
        tmp.status = TrackStatus::Confirmed;
        tmp.hits = track.hits;

        like[j] = measurement_likelihood(tmp, z, cfg);
        if (!track_update(tmp, z, cfg)) {
            like[j] = kZero;
            continue;
        }
        any = true;
        track.model_state[j] = tmp.state;
        for (std::size_t k = 0; k < 16; ++k) model_cov(track, j)[k] = tmp.covariance[k];
    }
    if (!any) return false;

    // mu_j proportional to c_j * Lambda_j.
    Real sum = kZero;
    for (std::size_t j = 0; j < kImmModels; ++j) {
        track.model_prob[j] = cbar[j] * like[j];
        sum += track.model_prob[j];
    }
    if (!(sum > kZero)) {
        // Every model found the measurement impossible. Falling back to the
        // prior is the only stable choice; zeroing would leave no filter at all.
        for (std::size_t j = 0; j < kImmModels; ++j) track.model_prob[j] = cbar[j];
        sum = kOne;
    }
    Real floored = kZero;
    for (std::size_t j = 0; j < kImmModels; ++j) {
        Real p = track.model_prob[j] / sum;
        if (p < imm.min_model_prob) p = imm.min_model_prob;
        track.model_prob[j] = p;
        floored += p;
    }
    for (std::size_t j = 0; j < kImmModels; ++j) track.model_prob[j] /= floored;

    imm_update_combined_estimate(track);
    track.last_update_s = z.time_s;
    track.label = z.label;
    track.amplitude = (track.hits == 0)
                    ? z.amplitude
                    : static_cast<Real>(0.7) * track.amplitude
                    + static_cast<Real>(0.3) * z.amplitude;
    ++track.hits;
    track.misses = 0;
    return true;
}

void imm_update_combined_estimate(ImmTrack& track) noexcept {
    Real xm[kN] = {kZero, kZero, kZero, kZero};
    for (std::size_t j = 0; j < kImmModels; ++j) {
        Real xj[kN];
        state_to_array(track.model_state[j], xj);
        for (std::size_t k = 0; k < kN; ++k) xm[k] += track.model_prob[j] * xj[k];
    }
    array_to_state(xm, track.state);

    for (std::size_t k = 0; k < 16; ++k) track.covariance[k] = kZero;
    for (std::size_t j = 0; j < kImmModels; ++j) {
        Real xj[kN];
        state_to_array(track.model_state[j], xj);
        Real d[kN];
        for (std::size_t k = 0; k < kN; ++k) d[k] = xj[k] - xm[k];
        const Real* pj = &track.model_cov[j * 16];
        for (std::size_t a = 0; a < kN; ++a) {
            for (std::size_t b = 0; b < kN; ++b) {
                track.covariance[a * kN + b] +=
                    track.model_prob[j] * (pj[a * kN + b] + d[a] * d[b]);
            }
        }
    }
}

Real imm_manoeuvre_probability(const ImmTrack& track) noexcept {
    return kOne - track.model_prob[0];
}

Real imm_turn_rate_estimate(const ImmTrack& track, const ImmConfig& imm) noexcept {
    Real w = kZero;
    for (std::size_t j = 0; j < kImmModels; ++j) {
        w += track.model_prob[j] * kModelRateScale[j] * imm.turn_rate_rps;
    }
    return w;
}



Real imm_track_nis(const ImmTrack& track, const Measurement& z,
                   const TrackerConfig& cfg) noexcept {
    if (!track.live()) return kHugeNis;
    Track tmp;
    tmp.state = track.state;
    tmp.covariance = track.covariance;
    tmp.status = TrackStatus::Confirmed;
    return track_nis(tmp, z, cfg);
}

std::size_t imm_count_tracks(std::span<const ImmTrack> tracks, TrackStatus status) noexcept {
    std::size_t n = 0;
    for (const ImmTrack& t : tracks) {
        if (t.status == status) ++n;
    }
    return n;
}

std::size_t imm_count_established(std::span<const ImmTrack> tracks) noexcept {
    std::size_t n = 0;
    for (const ImmTrack& t : tracks) {
        if (t.status == TrackStatus::Confirmed || t.status == TrackStatus::Coasting) ++n;
    }
    return n;
}

std::size_t imm_tracker_step(std::span<ImmTrack> tracks,
                             std::span<const Measurement> measurements,
                             const TrackerConfig& cfg,
                             const ImmConfig& imm,
                             Real time_s,
                             std::uint32_t& next_id) noexcept {
    // --- 1. Predict ---------------------------------------------------------
    for (ImmTrack& t : tracks) {
        if (!t.live()) continue;
        const Real dt = time_s - t.last_update_s;
        if (dt > kZero) imm_predict(t, dt, cfg, imm);
    }

    // --- 2. Association, in global cost order --------------------------------
    // Identical in structure to tracker_step's: build every gating pair, sort by
    // NIS, assign best-first. The cost is the COMBINED estimate's NIS, for the
    // reason given in the header -- a per-model gate lets the worst model veto.
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
            const Real d = imm_track_nis(tracks[ti], measurements[m], cfg);
            const Real gate = measurements[m].has_range_rate ? cfg.gate_chi2_3dof
                                                             : cfg.gate_chi2_2dof;
            if (d < gate) {
                pairs[n_pairs++] = Pair{d, static_cast<std::uint16_t>(ti),
                                        static_cast<std::uint16_t>(m)};
            }
        }
    }
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
        if (imm_update(tracks[ti], measurements[mi], cfg, imm)) {
            track_taken[ti] = true;
            used[mi] = true;
        }
    }
    for (std::size_t ti = 0; ti < tracks.size(); ++ti) {
        if (!tracks[ti].live()) continue;
        if (ti >= n_tracks || !track_taken[ti]) ++tracks[ti].misses;
    }

    // --- 3. Confirmation and deletion ----------------------------------------
    for (ImmTrack& t : tracks) {
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

    // --- 4. Initiate ---------------------------------------------------------
    for (std::size_t m = 0; m < n_meas; ++m) {
        if (used[m]) continue;
        for (ImmTrack& t : tracks) {
            if (t.live()) continue;
            imm_initiate(t, measurements[m], cfg, imm, next_id++);
            break;
        }
    }

    std::size_t live = 0;
    for (const ImmTrack& t : tracks) {
        if (t.live()) ++live;
    }
    return live;
}

}  // namespace phantom
