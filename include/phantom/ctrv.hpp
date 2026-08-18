// SPDX-License-Identifier: Apache-2.0
// libphantom-sonar — coordinated turn with the turn rate ESTIMATED.
//
// The IMM of v0.11 brackets a manoeuvre with models at 0 and +/-omega and
// reports a probability-weighted blend. That blend cannot leave [-omega, +omega]
// however hard the target turns, so a target turning at 6 deg/s against models
// set for 3 saturates the estimate and the filter under-predicts the turn.
//
// This filter estimates omega instead, by carrying it as a fifth state:
//
//     [ x, y, vx, vy, omega ]
//
// The price is that the transition stops being linear -- omega multiplies the
// velocity terms -- so the covariance must be propagated through a Jacobian
// rather than a constant matrix. That is the whole trade: an IMM is a bank of
// linear filters and a mixing rule, this is one nonlinear filter.
//
// WHICH TO USE. They are not ordered, and this one is not simply better:
//
//   IMM   converges faster after a manoeuvre STARTS, because a model that
//         already fits is waiting to take over. Cannot report a rate outside
//         its bracket. Robust: no model can diverge, only lose probability.
//
//   CTRV  reports the actual turn rate, with no ceiling. Slower: omega is
//         observed only through its effect on predicted position, so it takes
//         several scans of turning to be worth anything. Its uncertainty settles
//         to a steady state rather than shrinking -- see turn_rate_noise_rps,
//         where a filter that IS allowed to shrink it stops working.
//
// The honest summary is that CTRV measures a manoeuvre and an IMM reacts to
// one, and the test suite reports both rather than declaring a winner.
#ifndef PHANTOM_CTRV_HPP
#define PHANTOM_CTRV_HPP

#include "phantom/tracker.hpp"
#include "phantom/types.hpp"

#include <array>

namespace phantom {

inline constexpr std::size_t kCtrvN = 5;

struct CtrvState {
    Real x = 0;
    Real y = 0;
    Real vx = 0;
    Real vy = 0;
    Real turn_rate_rps = 0;

    [[nodiscard]] TargetState target_state() const noexcept {
        TargetState s;
        s.x = x; s.y = y; s.vx = vx; s.vy = vy;
        return s;
    }
};

struct CtrvConfig {
    // Unmodelled linear acceleration, as everywhere else in the library.
    Real process_accel_mps2 = static_cast<Real>(0.3);

    // Random walk on the turn rate itself, rad/s per sqrt(s). This is the knob
    // with no counterpart in a constant-velocity filter, and it is the one that
    // matters most.
    //
    // Setting it to ZERO is the trap, and it does not look like one: the
    // covariance shrinks, the filter grows confident, and the gain on omega
    // dies -- after which it cannot follow a change at all. Measured on a target
    // that starts turning at 5 deg/s, 30 scans later:
    //
    //     q_w      reported sigma    estimate (truth 5.00)
    //     0                0.11 deg/s           2.08      <- converged, wrong
    //     0.005            0.91                 4.97
    //     0.01             1.54                 4.90
    //     0.02             2.67                 4.89
    //
    // With q_w = 0 it converges beautifully to the wrong answer. The default
    // below is a compromise: 0.005 is more precise on a steady turn, but a
    // larger value notices a NEW manoeuvre sooner, and a tracker sees more new
    // manoeuvres than steady turns.
    Real turn_rate_noise_rps = static_cast<Real>(0.01);

    // Initial standard deviation of the turn rate. A single detection says
    // nothing about it, so this must cover the manoeuvres of interest.
    Real init_turn_rate_sigma_rps = deg2rad(static_cast<Real>(4));

    // Estimates are clamped here. Not cosmetic: the transition divides by
    // omega, and an unclamped filter that wanders to a large omega predicts a
    // tight circle, gates every measurement out and never recovers.
    Real max_turn_rate_rps = deg2rad(static_cast<Real>(20));
};

struct CtrvTrack {
    CtrvState state;
    std::array<Real, kCtrvN * kCtrvN> covariance{};

    Real last_update_s = 0;
    std::uint32_t id = 0;
    std::uint32_t hits = 0;
    std::uint32_t misses = 0;
    std::uint32_t age = 0;
    TrackStatus status = TrackStatus::Free;

    [[nodiscard]] bool live() const noexcept { return status != TrackStatus::Free; }
};

// Starts a track from one detection: position from the measurement, zero
// velocity and zero turn rate, with covariances wide enough that the next few
// updates are not gated out.
void ctrv_initiate(CtrvTrack& track, const Measurement& z,
                   const TrackerConfig& cfg, const CtrvConfig& ct,
                   std::uint32_t id) noexcept;

// Nonlinear prediction, with the covariance carried through the Jacobian.
//
// Both the transition and its Jacobian contain sin(wT)/w and (1-cos(wT))/w,
// which are 0/0 at w = 0. Below a threshold the implementation switches to the
// series expansions -- not as an approximation, but because they are the more
// accurate way to compute the same quantity there.
void ctrv_predict(CtrvTrack& track, Real dt,
                  const TrackerConfig& cfg, const CtrvConfig& ct) noexcept;

// EKF measurement update. Range and bearing (and range rate when present) do
// not depend on the turn rate, so the measurement Jacobian's fifth column is
// zero -- omega is observed only through its effect on the predicted position,
// which is why several scans of turning are needed before it is worth anything.
bool ctrv_update(CtrvTrack& track, const Measurement& z,
                 const TrackerConfig& cfg, const CtrvConfig& ct) noexcept;

// Normalised innovation squared, for gating.
[[nodiscard]] Real ctrv_nis(const CtrvTrack& track, const Measurement& z,
                            const TrackerConfig& cfg) noexcept;

// Standard deviation of the turn-rate estimate, rad/s. Read this before
// trusting the turn rate.
//
// It settles to a steady state set by turn_rate_noise_rps against the geometry
// -- measured at 1.5 deg/s for the default, on a 30 m/s target at 7 km with
// half-degree bearings. It does NOT shrink towards zero, and must not: omega
// can change at any moment. A filter reporting near-zero uncertainty about a
// turn rate has stopped being able to notice a new manoeuvre.
[[nodiscard]] Real ctrv_turn_rate_sigma_rps(const CtrvTrack& track) noexcept;

}  // namespace phantom

#endif  // PHANTOM_CTRV_HPP
