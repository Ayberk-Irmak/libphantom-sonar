// SPDX-License-Identifier: Apache-2.0
// libphantom-sonar — interacting multiple model filter.
//
// The tracker of v0.9 runs one constant-velocity model, and a CV filter has
// exactly two ways to handle a turn, both bad. Small process noise: the filter
// lags the manoeuvre, the innovations blow past the gate, and the track is
// dropped for a target that is plainly still there. Large process noise: the
// filter follows the turn, and spends the other 95% of the time chasing
// measurement noise with a velocity estimate too loose to be worth having.
//
// The IMM refuses the choice. It runs several models at once, keeps a
// probability for each, and mixes them -- so the quiet model does the work
// while the target flies straight, the turning models take over within a scan
// or two of a manoeuvre, and the estimate is the probability-weighted blend
// throughout. The model probabilities are also an output in their own right: a
// jump in the turning models' probability IS a manoeuvre detection, with no
// separate detector and no threshold to tune.
//
// Three models, because that is the smallest set that covers the problem:
//
//   0   CONSTANT VELOCITY, small process noise. The straight-and-level model.
//   1   COORDINATED TURN at +omega.
//   2   COORDINATED TURN at -omega.
//
// A coordinated turn at KNOWN rate is linear in the same four-state Cartesian
// vector the rest of the library uses, so all three models share a state space
// and no augmentation is needed. Estimating omega itself would need a fifth
// state and a nonlinear transition; two fixed rates bracket the manoeuvre
// instead, and the mixture handles rates in between. That is the standard
// trade, and it is why the turn-rate estimate this header exposes is a
// probability-weighted blend of the model rates rather than a measurement.
//
// SCOPE, stated plainly: this is the filter, verified on its own. It is not yet
// wired into tracker_step() -- association, gating and track management still
// run the single-model EKF. Joining them is v0.12.
#ifndef PHANTOM_IMM_HPP
#define PHANTOM_IMM_HPP

#include "phantom/tracker.hpp"
#include "phantom/types.hpp"

#include <array>

namespace phantom {

inline constexpr std::size_t kImmModels = 3;

struct ImmConfig {
    // Turn rate of the two manoeuvre models. 3 deg/s is a hard turn for a ship
    // and a gentle one for a torpedo; it should be set from the target class,
    // and setting it too small is the failure that looks like the IMM "not
    // working" -- the turning models then fit the manoeuvre no better than CV.
    Real turn_rate_rps = deg2rad(static_cast<Real>(3));

    // Process noise per model. The whole point is that these DIFFER: the CV
    // model is quiet so it estimates velocity well, the turning models are
    // noisy so they can absorb what the fixed turn rate does not capture.
    Real cv_accel_mps2 = static_cast<Real>(0.3);
    Real ct_accel_mps2 = static_cast<Real>(2.0);

    // Probability of leaving the current model in one scan, spread evenly over
    // the others. This sets how fast the mixture can switch: too small and the
    // IMM reacts as slowly as a single model, too large and it flickers.
    Real switch_probability = static_cast<Real>(0.05);

    // Probabilities are floored here after every update. A model that reaches
    // exactly zero can never come back, however well it would fit later.
    Real min_model_prob = static_cast<Real>(1e-4);
};

struct ImmTrack {
    // The combined estimate: the probability-weighted mixture, and the
    // covariance that accounts for the spread BETWEEN models as well as within
    // them. This is what a consumer should read.
    TargetState state;
    std::array<Real, 16> covariance{};

    // Per-model state, covariance (row-major 4x4 each) and probability.
    std::array<TargetState, kImmModels> model_state{};
    std::array<Real, kImmModels * 16> model_cov{};
    std::array<Real, kImmModels> model_prob{};

    Real last_update_s = 0;
    std::uint32_t hits = 0;
    bool live = false;
};

// Starts an IMM track from one detection. All models get the same state; the
// probabilities start on the CV model, since a target is presumed not to be
// manoeuvring until it shows otherwise.
void imm_initiate(ImmTrack& track, const Measurement& z,
                  const TrackerConfig& cfg, const ImmConfig& imm) noexcept;

// Mixing followed by per-model prediction.
//
// The mixing is what makes an IMM more than a bank of independent filters:
// before each model predicts, its state is replaced by a blend of ALL models'
// states, weighted by the probability that the target was in model i and
// switched to model j. A model that has been idle therefore starts from
// somewhere sensible when it takes over, instead of from wherever it drifted.
void imm_predict(ImmTrack& track, Real dt,
                 const TrackerConfig& cfg, const ImmConfig& imm) noexcept;

// Updates every model with the measurement, re-weights the models by how well
// each predicted it, and recombines. Returns false if no model could be
// updated, which means every one of them had a degenerate covariance.
bool imm_update(ImmTrack& track, const Measurement& z,
                const TrackerConfig& cfg, const ImmConfig& imm) noexcept;

// Recomputes `state` and `covariance` from the per-model estimates. Called
// automatically by imm_predict and imm_update; exposed because a caller that
// edits model probabilities by hand must restore the invariant.
void imm_update_combined_estimate(ImmTrack& track) noexcept;

// Probability that the target is turning: 1 - P(constant velocity).
//
// This is a manoeuvre detector that costs nothing extra, because the filter had
// to compute it anyway to do the mixing.
[[nodiscard]] Real imm_manoeuvre_probability(const ImmTrack& track) noexcept;

// Probability-weighted turn rate, rad/s. Signed: positive is the +omega model.
//
// NOT a measurement of the turn rate. With models at 0 and +/-omega it can only
// return a value in [-omega, +omega], and a target turning faster saturates it.
// Read it as "which way and roughly how hard", not as a rate.
[[nodiscard]] Real imm_turn_rate_estimate(const ImmTrack& track,
                                          const ImmConfig& imm) noexcept;

}  // namespace phantom

#endif  // PHANTOM_IMM_HPP
