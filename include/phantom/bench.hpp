// SPDX-License-Identifier: Apache-2.0
// libphantom-sonar — turning a recording into a measurement.
//
// v0.15 built an air bench and could not run it: the development machine has no
// working acoustic path. That failure is the reason this header exists, because
// the way it was found matters more than the fact of it.
//
// It was NOT found by looking at a recording and judging it. The recording
// looked plausible -- 13000 distinct sample values, a healthy RMS, a matched
// filter peak 44 dB above the background. Every one of those is consistent with
// a working microphone, and every one of them was produced by electrical noise
// and a startup transient. What settled it was a CONTROLLED comparison:
// recording while playing, recording in silence, and differencing the in-band
// energy. The answer came out NEGATIVE -- playing made it quieter -- which no
// working channel can do.
//
// So this header provides three things, in the order a real measurement needs
// them:
//
//   QUALIFY    decide whether a channel is carrying the signal at all, by the
//              A/B comparison above. Before this passes, no number taken from
//              that channel means anything, and the failure mode is a plausible
//              number rather than an obvious error.
//
//   LATENCY    remove the sound card's own delay. It is 10-50 ms where the
//              acoustic flight time over a desk is 1.5 ms, so it dominates
//              completely. The difference of two measurements at known
//              distances cancels it exactly.
//
//   RESPONSE   estimate the channel impulse response by regularised
//              deconvolution, so multipath can be compared against the model in
//              §10 rather than admired.
#ifndef PHANTOM_BENCH_HPP
#define PHANTOM_BENCH_HPP

#include "phantom/fft.hpp"
#include "phantom/types.hpp"

#include <span>

namespace phantom::bench {

// ---------------------------------------------------------------------------
// Channel qualification
// ---------------------------------------------------------------------------

struct QualifyResult {
    // In-band energy with the probe playing, and with silence, in dB relative
    // to an arbitrary common reference. Only the difference is meaningful.
    Real active_db = 0;
    Real silent_db = 0;
    // active_db - silent_db. A working channel makes this clearly positive.
    Real excess_db = 0;
    // Fraction of samples at or beyond full scale in the active recording.
    Real clipped_fraction = 0;
    // True only if the channel is carrying the probe AND is not clipping.
    bool usable = false;
};

struct QualifyConfig {
    Real sample_rate_hz = 48000;
    // The band the probe occupies. Energy outside it is not evidence.
    Real band_low_hz  = 2000;
    Real band_high_hz = 8000;
    // How much in-band excess counts as "the channel carries the probe".
    //
    // 6 dB is deliberately modest: the point is to catch a channel that is
    // carrying NOTHING, not to grade a good one. The dead channel this was
    // written against measured -4.60 dB, and no amount of gain moved it
    // positive.
    Real min_excess_db = 6;
    // Above this, the recording is clipping and its levels are fiction even if
    // the excess looks fine.
    Real max_clipped_fraction = static_cast<Real>(0.001);
};

// Compares a recording made while the probe played against one made in silence.
//
// `work` must be at least 2 * fft.size() complex. Returns false on bad input;
// on success the verdict is in `out`.
//
// BOTH recordings must come from the same device at the same gain, made close
// together in time. The whole method rests on the only difference being the
// probe, so a gain change between them invalidates it silently.
bool qualify_channel(const FftView& fft,
                     std::span<const Real> active,
                     std::span<const Real> silent,
                     const QualifyConfig& cfg,
                     std::span<Complex> work,
                     QualifyResult& out) noexcept;

// ---------------------------------------------------------------------------
// Latency cancellation
// ---------------------------------------------------------------------------

struct TwoDistanceResult {
    Real sound_speed_mps = 0;      // recovered from the two measurements
    Real fixed_latency_s = 0;      // the electronic delay, now separated out
    bool valid = false;
};

// Recovers both the sound speed and the fixed electronic latency from two
// delay measurements at two known distances.
//
//   t1 = d1/c + L
//   t2 = d2/c + L      =>   c = (d2 - d1)/(t2 - t1),  L = t1 - d1/c
//
// The subtraction is the whole point: L is 10-50 ms of sound-card buffering
// against 1.5 ms of flight time over half a metre, so a single measurement is
// 97% latency. Two measurements make the acoustics the only thing that differs.
//
// The distances must differ by enough that (t2 - t1) is well above the timing
// resolution -- returns valid = false otherwise rather than dividing by a
// number that is mostly noise.
[[nodiscard]] TwoDistanceResult two_distance_solve(Real distance1_m, Real delay1_s,
                                                   Real distance2_m, Real delay2_s,
                                                   Real timing_resolution_s) noexcept;

// ---------------------------------------------------------------------------
// Impulse response
// ---------------------------------------------------------------------------

// Estimates the channel impulse response by regularised deconvolution:
//
//   H = conj(S) R / (|S|^2 + eps * mean|S|^2)
//
// The regularisation is not optional. A probe has little energy outside its own
// band, so dividing by |S|^2 there amplifies whatever noise is present by an
// unbounded factor -- and the result still looks like an impulse response.
// `epsilon` sets the floor as a fraction of the probe's mean band power; 1e-3 is
// a reasonable start and larger values trade resolution for stability.
//
// `probe` and `recording` are both time domain and are zero-padded to the
// transform length. `work` must be at least 2 * fft.size() complex.
//
// Returns the number of impulse-response samples written, 0 on bad input.
std::size_t estimate_impulse_response(const FftView& fft,
                                      std::span<const Real> probe,
                                      std::span<const Real> recording,
                                      Real epsilon,
                                      std::span<Complex> work,
                                      std::span<Real> out) noexcept;

struct Arrival {
    Real delay_s = 0;
    Real amplitude = 0;
    Real relative_db = 0;   // relative to the strongest arrival
};

// Picks arrivals out of an impulse response.
//
// `min_separation_s` must be a few times the probe's own resolution (1/bandwidth)
// and no more. An earlier version of the bench tool used 20 ms, which on a desk
// discarded every echo -- the direct path and its first reflection are typically
// 2-4 ms apart, and 20 ms is seven metres of extra path.
std::size_t find_arrivals(std::span<const Real> impulse_response,
                          Real sample_rate_hz,
                          Real min_separation_s,
                          Real threshold_db,
                          std::span<Arrival> out) noexcept;

}  // namespace phantom::bench

#endif  // PHANTOM_BENCH_HPP
