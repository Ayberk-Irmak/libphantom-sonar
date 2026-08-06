// SPDX-License-Identifier: Apache-2.0
#include "phantom/waveform.hpp"

#include <cmath>

namespace phantom {
namespace {

constexpr Real kZero = static_cast<Real>(0);
constexpr Real kOne  = static_cast<Real>(1);
constexpr Real kTwo  = static_cast<Real>(2);

// Below this fractional sweep an HFM is indistinguishable from an LFM in
// double precision, and its closed-form phase loses significance to
// cancellation. The linear form is used instead; the two agree to O(k^2).
constexpr Real kHfmLinearEps = static_cast<Real>(1e-9);

}  // namespace

const char* pulse_type_name(PulseType t) noexcept {
    switch (t) {
        case PulseType::Cw:      return "CW";
        case PulseType::LfmUp:   return "LFM-up";
        case PulseType::LfmDown: return "LFM-down";
        case PulseType::Hfm:     return "HFM";
        case PulseType::Unknown: break;
    }
    return "unknown";
}

Real PulseSpec::centre_frequency_hz() const noexcept {
    if (type == PulseType::Cw) return f_start_hz;
    return (f_start_hz + f_end_hz) / kTwo;
}

Real PulseSpec::bandwidth_hz() const noexcept {
    if (type == PulseType::Cw) return kZero;
    return std::fabs(f_end_hz - f_start_hz);
}

Real PulseSpec::chirp_rate_hz_s() const noexcept {
    if (type == PulseType::Cw || !(duration_s > kZero)) return kZero;
    return (f_end_hz - f_start_hz) / duration_s;
}

Real PulseSpec::time_bandwidth_product() const noexcept {
    return bandwidth_hz() * duration_s;
}

bool PulseSpec::valid() const noexcept {
    if (!(duration_s > kZero) || !(amplitude > kZero)) return false;
    if (!(f_start_hz > kZero)) return false;
    switch (type) {
        case PulseType::Cw:
            return true;
        case PulseType::LfmUp:
            return f_end_hz > f_start_hz;
        case PulseType::LfmDown:
            return f_end_hz < f_start_hz && f_end_hz > kZero;
        case PulseType::Hfm:
            // An HFM's phase is built from 1/f, so neither endpoint may be zero
            // and they must differ for the hyperbolic form to mean anything.
            return f_end_hz > kZero && std::fabs(f_end_hz - f_start_hz) > kZero;
        case PulseType::Unknown:
            break;
    }
    return false;
}

// phi(t) = 2*pi * integral of f(t) dt, with phi(0) = 0.
Real pulse_phase(const PulseSpec& spec, Real t) noexcept {
    const Real two_pi = kTwo * kPi;
    switch (spec.type) {
        case PulseType::Cw:
            return two_pi * spec.f_start_hz * t;

        case PulseType::LfmUp:
        case PulseType::LfmDown: {
            // f(t) = f0 + mu*t   ->   phi = 2*pi*(f0*t + mu*t^2/2)
            const Real mu = spec.chirp_rate_hz_s();
            return two_pi * (spec.f_start_hz * t + static_cast<Real>(0.5) * mu * t * t);
        }

        case PulseType::Hfm: {
            // 1/f(t) = 1/f0 + k*t   with   k = (1/f1 - 1/f0) / T
            // phi = 2*pi * integral dt / (1/f0 + k t) = (2*pi/k) * ln(1 + k*f0*t)
            const Real inv0 = kOne / spec.f_start_hz;
            const Real inv1 = kOne / spec.f_end_hz;
            const Real k = (inv1 - inv0) / spec.duration_s;
            if (std::fabs(k * spec.f_start_hz * spec.duration_s) < kHfmLinearEps) {
                return two_pi * spec.f_start_hz * t;
            }
            return (two_pi / k) * std::log(kOne + k * spec.f_start_hz * t);
        }

        case PulseType::Unknown:
            break;
    }
    return kZero;
}

Real pulse_frequency(const PulseSpec& spec, Real t) noexcept {
    switch (spec.type) {
        case PulseType::Cw:
            return spec.f_start_hz;
        case PulseType::LfmUp:
        case PulseType::LfmDown:
            return spec.f_start_hz + spec.chirp_rate_hz_s() * t;
        case PulseType::Hfm: {
            const Real inv0 = kOne / spec.f_start_hz;
            const Real inv1 = kOne / spec.f_end_hz;
            const Real k = (inv1 - inv0) / spec.duration_s;
            const Real inv = inv0 + k * t;
            return (inv > kZero) ? kOne / inv : kZero;
        }
        case PulseType::Unknown:
            break;
    }
    return kZero;
}

Real taper_weight(Taper taper, Real t, Real duration_s) noexcept {
    if (!(duration_s > kZero)) return kZero;
    const Real x = t / duration_s;  // normalised position in [0, 1]
    if (x < kZero || x > kOne) return kZero;

    switch (taper) {
        case Taper::Rectangular:
            return kOne;
        case Taper::Hann:
            return static_cast<Real>(0.5) * (kOne - std::cos(kTwo * kPi * x));
        case Taper::Tukey25: {
            // Cosine-tapered edges over the outer 12.5% at each end.
            constexpr Real alpha = static_cast<Real>(0.25);
            const Real half = alpha / kTwo;
            if (x < half) {
                return static_cast<Real>(0.5) * (kOne - std::cos(kPi * x / half));
            }
            if (x > kOne - half) {
                return static_cast<Real>(0.5) * (kOne - std::cos(kPi * (kOne - x) / half));
            }
            return kOne;
        }
    }
    return kOne;
}

std::size_t pulse_length(const PulseSpec& spec, Real sample_rate_hz) noexcept {
    if (!spec.valid() || !(sample_rate_hz > kZero)) return 0;
    const Real n = spec.duration_s * sample_rate_hz;
    if (!(n >= kOne)) return 0;
    return static_cast<std::size_t>(n);
}

std::size_t render_real(const PulseSpec& spec, Real sample_rate_hz,
                        std::span<Real> out) noexcept {
    const std::size_t n = pulse_length(spec, sample_rate_hz);
    if (n == 0 || out.size() < n) return 0;
    const Real dt = kOne / sample_rate_hz;
    for (std::size_t i = 0; i < n; ++i) {
        const Real t = static_cast<Real>(i) * dt;
        out[i] = spec.amplitude * taper_weight(spec.taper, t, spec.duration_s)
               * std::cos(pulse_phase(spec, t));
    }
    return n;
}

std::size_t render_analytic(const PulseSpec& spec, Real sample_rate_hz,
                            std::span<Complex> out) noexcept {
    const std::size_t n = pulse_length(spec, sample_rate_hz);
    if (n == 0 || out.size() < n) return 0;
    const Real dt = kOne / sample_rate_hz;
    for (std::size_t i = 0; i < n; ++i) {
        const Real t = static_cast<Real>(i) * dt;
        const Real a = spec.amplitude * taper_weight(spec.taper, t, spec.duration_s);
        const Real ph = pulse_phase(spec, t);
        out[i] = Complex(a * std::cos(ph), a * std::sin(ph));
    }
    return n;
}

std::size_t render_real_doppler(const PulseSpec& spec, Real sample_rate_hz,
                                Real doppler, std::span<Real> out) noexcept {
    if (!spec.valid() || !(sample_rate_hz > kZero)) return 0;
    const Real scale = kOne + doppler;
    if (!(scale > kZero)) return 0;

    // A closing target compresses the pulse in time; the received duration is
    // T/scale and the received waveform is s(scale * t).
    const Real rx_duration = spec.duration_s / scale;
    const Real n_real = rx_duration * sample_rate_hz;
    if (!(n_real >= kOne)) return 0;
    const auto n = static_cast<std::size_t>(n_real);
    if (out.size() < n) return 0;

    const Real dt = kOne / sample_rate_hz;
    for (std::size_t i = 0; i < n; ++i) {
        const Real t = static_cast<Real>(i) * dt;
        const Real t_src = scale * t;  // time in the transmitted waveform
        out[i] = spec.amplitude * taper_weight(spec.taper, t_src, spec.duration_s)
               * std::cos(pulse_phase(spec, t_src));
    }
    return n;
}

}  // namespace phantom
