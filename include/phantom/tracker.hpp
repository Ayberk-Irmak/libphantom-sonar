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
// A NOTE ON WHAT TRACKING DOES AND DOES NOT FIX. The roadmap for this release
// claimed that time consistency would suppress the cross-template ghosts of
// v0.2. That claim was wrong, and the test suite demonstrates it rather than
// quietly dropping it: a ghost appears whenever the real arrival does, at a
// fixed offset, so it is exactly as consistent over time as the target and
// forms its own perfectly healthy track. What tracking DOES kill is false
// alarms, which do not repeat -- and it kills them decisively.
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
    // Carried through so a track can report what kind of thing it is tracking;
    // not used by the filter.
    std::uint16_t label = 0;
};

struct TrackerConfig {
    // Standard deviation of the unmodelled acceleration, m/s^2. This is the
    // single knob that sets how much manoeuvre the filter expects: too small
    // and it lags a turn, too large and it chases noise.
    Real process_accel_mps2 = static_cast<Real>(0.5);

    Real range_sigma_m = 5;
    Real bearing_sigma_rad = deg2rad(static_cast<Real>(1));

    // Chi-square gate on the normalised innovation squared, 2 degrees of
    // freedom. 9.21 admits 99% of true measurements; 5.99 admits 95%.
    Real gate_chi2 = static_cast<Real>(9.21);

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
// d^2 = y^T S^-1 y, which is chi-square with 2 degrees of freedom when the
// filter is consistent. Returns a large value rather than a NaN for a
// degenerate track.
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

// The 95th and 99th percentiles of the chi-square distribution with 2 degrees
// of freedom, for choosing a gate: chi2 = -2 ln(1 - p).
[[nodiscard]] Real chi2_gate_2dof(Real probability) noexcept;

}  // namespace phantom

#endif  // PHANTOM_TRACKER_HPP
