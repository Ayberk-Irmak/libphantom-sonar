// SPDX-License-Identifier: Apache-2.0
// libphantom-sonar — active sonar pulse synthesis.
//
// Three pulse families, because they are what active sonars actually transmit
// and they behave very differently at the receiver:
//
//   CW   constant frequency. Excellent Doppler resolution, poor range
//        resolution (1/T), trivially detected.
//   LFM  frequency linear in time. Range resolution 1/B independent of pulse
//        length, at the cost of range-Doppler coupling.
//   HFM  *period* linear in time. Under a Doppler time-scaling an HFM maps
//        onto a delayed copy of itself, so the matched filter peak survives
//        rather than smearing. This is why HFM dominates underwater, where
//        v/c is ~500x larger than in radar: 1 m/s of closing speed is already
//        6.7e-4 of the wave speed, and it scales the waveform in time rather
//        than merely shifting its frequency.
//
// Phase is evaluated in closed form at each sample rather than accumulated,
// so there is no phase drift to integrate over a long pulse.
#ifndef PHANTOM_WAVEFORM_HPP
#define PHANTOM_WAVEFORM_HPP

#include "phantom/fft.hpp"
#include "phantom/types.hpp"

#include <span>

namespace phantom {

enum class PulseType : std::uint8_t {
    Unknown = 0,
    Cw,       // continuous wave
    LfmUp,    // linear FM, upsweep
    LfmDown,  // linear FM, downsweep
    Hfm,      // hyperbolic FM (upsweep when f_end > f_start)
};

[[nodiscard]] const char* pulse_type_name(PulseType t) noexcept;

// Amplitude taper. A rectangular pulse has -13 dB matched-filter sidelobes,
// which is enough to hide a weak second target beside a strong one; Hann drops
// them to about -32 dB at the cost of a wider mainlobe and ~1.8 dB of SNR.
enum class Taper : std::uint8_t { Rectangular = 0, Hann, Tukey25 };

struct PulseSpec {
    PulseType type = PulseType::Cw;
    Real f_start_hz = 0;
    Real f_end_hz = 0;      // ignored for Cw
    Real duration_s = 0;
    Real amplitude = 1;
    Taper taper = Taper::Rectangular;

    [[nodiscard]] Real centre_frequency_hz() const noexcept;
    [[nodiscard]] Real bandwidth_hz() const noexcept;
    // Chirp rate in Hz/s. Zero for CW; for HFM this is the mean rate, since an
    // HFM's instantaneous rate varies across the sweep by construction.
    [[nodiscard]] Real chirp_rate_hz_s() const noexcept;
    [[nodiscard]] Real time_bandwidth_product() const noexcept;
    [[nodiscard]] bool valid() const noexcept;
};

// Instantaneous phase in radians at time `t` seconds from the pulse start.
[[nodiscard]] Real pulse_phase(const PulseSpec& spec, Real t) noexcept;

// Instantaneous frequency in Hz at time `t`.
[[nodiscard]] Real pulse_frequency(const PulseSpec& spec, Real t) noexcept;

// Amplitude taper weight at time `t`, in [0, 1].
[[nodiscard]] Real taper_weight(Taper taper, Real t, Real duration_s) noexcept;

// Number of samples a pulse occupies at `sample_rate_hz`.
[[nodiscard]] std::size_t pulse_length(const PulseSpec& spec, Real sample_rate_hz) noexcept;

// Renders the real passband waveform into `out`, starting at sample 0.
// Returns the number of samples written (0 on invalid input).
std::size_t render_real(const PulseSpec& spec, Real sample_rate_hz,
                        std::span<Real> out) noexcept;

// Renders the analytic (complex) waveform exp(j*phase), used as the matched
// filter replica. Correlating a real signal against this yields the complex
// envelope directly, with no separate basebanding stage.
std::size_t render_analytic(const PulseSpec& spec, Real sample_rate_hz,
                            std::span<Complex> out) noexcept;

// Renders the waveform with its time axis scaled by `1 + doppler`, i.e. the
// waveform a receiver sees from a target closing at v = doppler * c.
//
// This is a time *dilation*, not a frequency shift. Modelling it as a shift is
// the standard mistake carried over from radar, and it is wrong underwater:
// with c ~ 1500 m/s the scaling is large enough to destroy an LFM's matched
// filter peak outright.
std::size_t render_real_doppler(const PulseSpec& spec, Real sample_rate_hz,
                                Real doppler, std::span<Real> out) noexcept;

// The analytic form of the same, used as a Doppler-bin replica.
std::size_t render_analytic_doppler(const PulseSpec& spec, Real sample_rate_hz,
                                    Real doppler, std::span<Complex> out) noexcept;

// The Doppler mismatch |delta| = |v|/c at which the matched filter peak drops
// by `max_loss_db`. This is what sets how many Doppler bins a bank needs.
//
// Derived per family, because they behave completely differently:
//
//   CW   the output is sinc(delta * f0 * T); expanding to second order,
//        delta = sqrt(6 (1 - a)) / (pi f0 T),   a = 10^(-L/20)
//
//   LFM  a mismatched chirp leaves a residual quadratic phase B t^2 with
//        B = mu delta. Removing its least-squares linear part (the delay the
//        detector picks up as a range bias) leaves an RMS phase of
//        2 pi TB delta / sqrt(180), and loss_dB = 4.343 sigma^2. Hence
//        delta_1dB is almost exactly 1/TB.
//
//   HFM  time scaling maps the family onto itself, so there is no phase
//        residual at all -- only the duration mismatch, an amplitude factor of
//        (1 - |delta|). Two orders of magnitude more tolerant than an LFM of
//        the same time-bandwidth product, which is the whole reason it is used.
//
// VALIDITY: these are small-mismatch expansions. They are accurate to a few
// dB of loss -- the regime a bank is actually designed around -- and are NOT a
// model of deep mismatch. Beyond ~5 dB the true loss saturates while the LFM
// expression runs away (it predicts 49 dB where 9.9 dB is measured). Use it to
// space bins, not to predict what happens outside the bank.
[[nodiscard]] Real doppler_tolerance(const PulseSpec& spec, Real max_loss_db) noexcept;

}  // namespace phantom

#endif  // PHANTOM_WAVEFORM_HPP
