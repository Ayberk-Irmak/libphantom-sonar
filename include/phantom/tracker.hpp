// SPDX-License-Identifier: Apache-2.0
// libphantom-sonar — detection-to-track association and target state.
//
// A Pulse Descriptor Word now carries time, type, Doppler and (since v0.8)
// bearing. Nothing has yet connected them across blocks, which is what turns a
// list of detections into a picture.
//
// The filter is an EKF over a constant-velocity target in Cartesian
// coordinates, measuring range and bearing. Cartesian because the dynamics are
// then linear; polar measurements because that is what a sonar produces. The
// nonlinearity lives entirely in the measurement Jacobian.
//
// A NOTE ON WHAT TRACKING DOES AND DOES NOT FIX. v0.8's roadmap claimed that
// time consistency would suppress the cross-template ghosts of v0.2. That was
// wrong: a ghost appears whenever the real arrival does, at a fixed offset, so
// it is exactly as consistent over time as the target and forms its own
// perfectly healthy track. Tracking kills FALSE ALARMS, which do not repeat.
//
// v0.10 kills the ghosts too, but by recognising them as template artefacts
// rather than by any appeal to time -- see suppress_template_ghosts below.
#ifndef PHANTOM_TRACKER_HPP
#define PHANTOM_TRACKER_HPP

#include "phantom/types.hpp"

#include <array>
#include <span>

namespace phantom {

// Geometry convention: y is along broadside (the direction a zero-bearing
// target lies in), x is cross-range. So bearing = atan2(x, y), matching the
// array convention of measuring from broadside.
struct TargetState {
    Real x = 0;    // cross-range, m
    Real y = 0;    // along-broadside range, m
    Real vx = 0;   // m/s
    Real vy = 0;   // m/s

    [[nodiscard]] Real range_m() const noexcept;
    [[nodiscard]] Real bearing_rad() const noexcept;
    // Closing rate, positive when the range is shrinking -- the sign the
    // Doppler bank reports.
    [[nodiscard]] Real range_rate_mps() const noexcept;
};

enum class TrackStatus : std::uint8_t { Free = 0, Tentative, Confirmed, Coasting };

struct Track {
    TargetState state;
    // 4x4 covariance, row-major, in the state's own units.
    std::array<Real, 16> covariance{};

    Real last_update_s = 0;
    // Smoothed measurement amplitude and the waveform label most recently
    // associated. Not used by the filter; used by ghost recognition.
    Real amplitude = 0;
    std::uint16_t label = 0;
    std::uint32_t id = 0;
    std::uint32_t hits = 0;        // updates since initiation
    std::uint32_t misses = 0;      // consecutive predictions without an update
    std::uint32_t age = 0;         // total steps
    TrackStatus status = TrackStatus::Free;

    [[nodiscard]] bool live() const noexcept { return status != TrackStatus::Free; }
};

// One detection reduced to what a tracker needs.
struct Measurement {
    Real range_m = 0;
    Real bearing_rad = 0;
    Real time_s = 0;

    // Closing rate from the Doppler bank, positive when the range is shrinking.
    // Until v0.10 this was measured and thrown away: the filter inferred
    // velocity from position history alone while a direct measurement of its
    // radial component sat unused in the PulseDescriptor.
    //
    // Set `has_range_rate` only when the bank actually resolved it -- a
    // zero-Doppler template reports 0 m/s, which is a very different statement
    // from "the target is stationary", and feeding it in as a measurement would
    // pin every track's radial velocity to zero.
    Real range_rate_mps = 0;
    bool has_range_rate = false;

    // Which waveform matched. Ghost recognition needs it: a cross-template
    // ghost is by definition a DIFFERENT template from the arrival that caused
    // it, and two real targets lit by one sonar return the same waveform.
    std::uint16_t label = 0;
    // Detected amplitude, for the ghost's amplitude-ratio test.
    Real amplitude = 1;
};

struct TrackerConfig {
    // Standard deviation of the unmodelled acceleration, m/s^2. This is the
    // single knob that sets how much manoeuvre the filter expects: too small
    // and it lags a turn, too large and it chases noise.
    Real process_accel_mps2 = static_cast<Real>(0.5);

    Real range_sigma_m = 5;
    Real bearing_sigma_rad = deg2rad(static_cast<Real>(1));
    // Standard deviation of the Doppler bank's closing-rate estimate. Half a
    // bin spacing is the right order: the bank quantises, so this is not a
    // noise figure so much as a quantisation one.
    Real range_rate_sigma_mps = 3;

    // Chi-square gate on the normalised innovation squared. The dimension
    // depends on the measurement: 2 without a range rate, 3 with one, so a
    // fixed threshold means different gate probabilities for the two. Use
    // chi2_gate() with the matching dof.
    Real gate_chi2_2dof = static_cast<Real>(9.21);
    Real gate_chi2_3dof = static_cast<Real>(11.34);

    // A track confirms on `confirm_hits` updates and is dropped after
    // `delete_misses` consecutive misses.
    std::uint32_t confirm_hits = 3;
    std::uint32_t delete_misses = 3;

    // Initial velocity uncertainty for a new track, m/s. A single detection
    // says nothing about velocity, so this must be large enough to cover the
    // targets of interest or the first update will be gated out.
    Real init_velocity_sigma_mps = 10;
};

// Predicts a track forward by `dt` seconds under the constant-velocity model,
// growing the covariance by the discrete white-noise-acceleration Q.
void track_predict(Track& track, Real dt, const TrackerConfig& cfg) noexcept;

// Normalised innovation squared between a predicted track and a measurement:
// d^2 = y^T S^-1 y, chi-square with 2 degrees of freedom -- or 3 when the
// measurement carries a range rate. Returns a large value rather than a NaN
// for a degenerate track.
//
// This is also the gating statistic, and the quantity whose distribution the
// test suite checks -- a filter that tracks well but reports the wrong
// covariance passes every position test and fails this one.
[[nodiscard]] Real track_nis(const Track& track, const Measurement& z,
                             const TrackerConfig& cfg) noexcept;

// EKF measurement update. Returns false if the innovation covariance is
// singular, which means the track's own covariance had already collapsed.
bool track_update(Track& track, const Measurement& z, const TrackerConfig& cfg) noexcept;

// Starts a track from a single detection: position from the measurement, zero
// velocity, and a velocity covariance wide enough that the next update is not
// gated out.
void track_initiate(Track& track, const Measurement& z, const TrackerConfig& cfg,
                    std::uint32_t id) noexcept;

// One tracker step: predict every live track to `time_s`, associate each
// measurement to its nearest gating track, update, then apply the confirmation
// and deletion rules and start tracks from whatever was left over.
//
// Association is greedy nearest-neighbour on the NIS. That is not optimal --
// a global assignment would do better when two targets cross -- but it is
// O(tracks * measurements) with no allocation, and the failure mode (a swap
// during a crossing) is well understood rather than surprising.
//
// Returns the number of live tracks after the step.
std::size_t tracker_step(std::span<Track> tracks,
                         std::span<const Measurement> measurements,
                         const TrackerConfig& cfg,
                         Real time_s,
                         std::uint32_t& next_id) noexcept;

// Counts tracks in a given state.
[[nodiscard]] std::size_t count_tracks(std::span<const Track> tracks,
                                       TrackStatus status) noexcept;

// Confirmed plus Coasting: the tracks a caller would call real.
//
// Coasting means "confirmed, but missed the most recent scan", which is a
// statement about the last update and not about whether the target exists.
// Counting only Confirmed undercounts by exactly the tracks that happened to be
// missed this instant, which for a detector running at any realistic Pd is a
// large fraction of them.
[[nodiscard]] std::size_t count_established(std::span<const Track> tracks) noexcept;

// Gaussian likelihood of a measurement given a predicted track:
//
//   L = exp(-d^2/2) / sqrt( (2 pi)^m |S| )
//
// where d^2 is the NIS and S the innovation covariance. Returns 0 for a
// degenerate track. This is what an IMM weighs its models by -- the NIS alone
// is not enough, because a model can win on residual simply by being vaguer,
// and the |S| term is what charges it for that.
[[nodiscard]] Real measurement_likelihood(const Track& track, const Measurement& z,
                                          const TrackerConfig& cfg) noexcept;

// Chi-square quantile for choosing a gate.
//
// 2 dof inverts in closed form (the CDF is 1 - exp(-x/2), so x = -2 ln(1-p)).
// 3 dof does not, so it is bisected on the closed-form CDF
// erf(sqrt(x/2)) - sqrt(2x/pi) exp(-x/2).
[[nodiscard]] Real chi2_gate_2dof(Real probability) noexcept;
[[nodiscard]] Real chi2_gate(Real probability, std::size_t dof) noexcept;

// ---------------------------------------------------------------------------
// Cross-template ghost recognition
// ---------------------------------------------------------------------------
//
// v0.9 measured that tracking cannot suppress these: a ghost is exactly as
// consistent over time as the target that produces it. What distinguishes it is
// not its kinematics but its ORIGIN -- it is the same arrival seen through a
// different matched filter, so it shares the target's bearing and bearing rate,
// sits at a fixed range offset, is weaker, and carries a DIFFERENT waveform
// label.
//
// The label is what makes this safe. Two real targets illuminated by one sonar
// return the same waveform, so they share a label and are never paired. Without
// that check, two targets in line astern -- same bearing, same kinematics,
// fixed separation -- would be indistinguishable from a ghost pair, and the
// test suite demonstrates exactly that.
struct GhostConfig {
    // Ghosts are the same arrival, so they share a TRUE bearing exactly -- but
    // the two tracks estimate it independently, so the tolerance must cover a
    // few times the achieved bearing accuracy, not the geometry. Two tracks on
    // one true bearing with 1 degree measurements differ by 1.7 degrees often
    // enough to matter; set this from what the array actually delivers.
    //
    // Loosening it is safe because the LABEL check, not this one, is what
    // separates a ghost from a real contact.
    Real max_bearing_delta_rad = deg2rad(static_cast<Real>(3));
    // ... and move together.
    Real max_range_rate_delta_mps = 2;
    // The offset must be a plausible template cross-correlation lag, not a
    // second target at some arbitrary range.
    Real max_range_offset_m = 400;
    // The ghost must be the weaker of the pair by at least this much.
    Real min_amplitude_ratio_db = 4;
    // Both tracks must have been seen this often, so a transient pairing does
    // not delete a real track.
    std::uint32_t min_hits = 4;
};

// Marks the weaker member of each recognised ghost pair as Free. Returns the
// number suppressed.
//
// Only pairs with DIFFERENT labels are considered, which is the whole safety
// argument -- see the note above.
std::size_t suppress_template_ghosts(std::span<Track> tracks,
                                     const GhostConfig& cfg) noexcept;

}  // namespace phantom

#endif  // PHANTOM_TRACKER_HPP
