// SPDX-License-Identifier: Apache-2.0
#include "phantom/array.hpp"

#include <cmath>

namespace phantom {
namespace {

constexpr Real kZero = static_cast<Real>(0);
constexpr Real kOne  = static_cast<Real>(1);
constexpr Real kTwo  = static_cast<Real>(2);

// Below this the array factor's 0/0 limit is taken as unity rather than
// evaluated. sin(psi/2) at 1e-9 still has ten digits of headroom in double and
// four in float, so the branch never fires on anything a caller can resolve.
constexpr Real kPsiEps = static_cast<Real>(1e-9);

}  // namespace

Real array_factor(const LineArray& array, Real lambda_m,
                  Real steer_rad, Real look_rad) noexcept {
    if (!array.valid() || !(lambda_m > kZero)) return kZero;
    const Real k = kTwo * kPi / lambda_m;
    const Real psi = k * array.spacing_m * (std::sin(look_rad) - std::sin(steer_rad));
    const Real half = psi / kTwo;
    const Real den = std::sin(half);
    const auto n = static_cast<Real>(array.element_count);
    if (std::fabs(den) < kPsiEps) return kOne;   // the 0/0 peak
    return std::fabs(std::sin(n * half) / (n * den));
}

Real first_null_offset_rad(const LineArray& array, Real lambda_m, Real steer_rad) noexcept {
    if (!array.valid() || !(lambda_m > kZero)) return kZero;
    const Real n = static_cast<Real>(array.element_count);
    // psi = 2 pi / N  ->  sin(look) - sin(steer) = lambda / (N d)
    const Real delta_sin = lambda_m / (n * array.spacing_m);
    const Real target = std::sin(steer_rad) + delta_sin;
    if (!(std::fabs(target) <= kOne)) return kZero;   // null is off the visible region
    return std::asin(target) - steer_rad;
}

Real beamwidth_3db_rad(const LineArray& array, Real lambda_m, Real steer_rad) noexcept {
    if (!array.valid() || !(lambda_m > kZero)) return kZero;
    const Real cos_th = std::fabs(std::cos(steer_rad));
    if (!(cos_th > static_cast<Real>(1e-6))) return kPi;   // endfire: the beam is a cone
    const Real n = static_cast<Real>(array.element_count);
    return static_cast<Real>(0.886) * lambda_m / (n * array.spacing_m * cos_th);
}

Real max_spacing_no_grating_lobes_m(Real lambda_m, Real max_steer_rad) noexcept {
    if (!(lambda_m > kZero)) return kZero;
    return lambda_m / (kOne + std::fabs(std::sin(max_steer_rad)));
}

bool has_grating_lobe(const LineArray& array, Real lambda_m, Real steer_rad) noexcept {
    if (!array.valid() || !(lambda_m > kZero)) return false;
    // A grating lobe appears where psi = +/-2 pi, i.e.
    // sin(look) = sin(steer) -/+ lambda/d, and it is real only if that lands
    // inside [-1, 1].
    const Real shift = lambda_m / array.spacing_m;
    const Real s = std::sin(steer_rad);
    return (std::fabs(s + shift) <= kOne) || (std::fabs(s - shift) <= kOne);
}

Real array_gain_db(const LineArray& array) noexcept {
    if (!array.valid()) return kZero;
    return static_cast<Real>(10) * std::log10(static_cast<Real>(array.element_count));
}

Real narrowband_bandwidth_limit_hz(const LineArray& array, Real steer_rad,
                                   Real sound_speed_mps) noexcept {
    if (!array.valid() || !(sound_speed_mps > kZero)) return kZero;
    const Real traversal = array.aperture_m() * std::fabs(std::sin(steer_rad)) / sound_speed_mps;
    if (!(traversal > kZero)) return kZero;   // broadside: no traversal, no limit
    return kOne / (static_cast<Real>(10) * traversal);
}

Real bearing_crlb_rad(const LineArray& array, Real lambda_m,
                      Real bearing_rad, Real element_snr) noexcept {
    if (!array.valid() || !(lambda_m > kZero) || !(element_snr > kZero)) return kZero;
    const Real cos_th = std::fabs(std::cos(bearing_rad));
    // 1e-6 rather than something tighter: cos(pi/2) evaluates to 6e-17 in
    // double but -4e-8 in float, so a 1e-9 guard silently stops firing in a
    // single-precision build and the bound comes back finite at endfire, where
    // there is no bearing information at all. Same precedent as
    // spreading_loss_db in v0.4.
    if (!(cos_th > static_cast<Real>(1e-6))) return kZero;   // endfire: no information

    const Real k = kTwo * kPi / lambda_m;
    const auto n = static_cast<Real>(array.element_count);
    // Sum of centred element indices squared, N(N^2-1)/12.
    const Real spread = n * (n * n - kOne) / static_cast<Real>(12);
    const Real kd_cos = k * array.spacing_m * cos_th;
    const Real fisher = kTwo * element_snr * kd_cos * kd_cos * spread;
    if (!(fisher > kZero)) return kZero;
    return kOne / std::sqrt(fisher);
}

std::size_t synthesize_plane_wave(const LineArray& array, Real lambda_m,
                                  Real bearing_rad, Real amplitude,
                                  std::span<Complex> out) noexcept {
    if (!array.valid() || !(lambda_m > kZero)) return 0;
    if (out.size() < array.element_count) return 0;

    const Real k = kTwo * kPi / lambda_m;
    const Real s = std::sin(bearing_rad);
    // Positions centred on the array, so the phase reference is its middle.
    const Real centre = static_cast<Real>(array.element_count - 1) / kTwo;
    for (std::size_t i = 0; i < array.element_count; ++i) {
        const Real x = (static_cast<Real>(i) - centre) * array.spacing_m;
        const Real ph = k * x * s;
        out[i] = Complex(amplitude * std::cos(ph), amplitude * std::sin(ph));
    }
    return array.element_count;
}

std::size_t beamform_power(const LineArray& array, Real lambda_m,
                           std::span<const Complex> elements,
                           Real angle_min_rad, Real angle_max_rad,
                           std::span<Real> out_power) noexcept {
    if (!array.valid() || !(lambda_m > kZero)) return 0;
    if (elements.size() < array.element_count || out_power.empty()) return 0;
    if (!(angle_max_rad > angle_min_rad)) return 0;

    const Real k = kTwo * kPi / lambda_m;
    const Real centre = static_cast<Real>(array.element_count - 1) / kTwo;
    const std::size_t n_angles = out_power.size();
    const Real step = (n_angles > 1)
                    ? (angle_max_rad - angle_min_rad) / static_cast<Real>(n_angles - 1)
                    : kZero;

    for (std::size_t a = 0; a < n_angles; ++a) {
        const Real steer = angle_min_rad + step * static_cast<Real>(a);
        const Real s = std::sin(steer);
        Real re = kZero;
        Real im = kZero;
        for (std::size_t i = 0; i < array.element_count; ++i) {
            const Real x = (static_cast<Real>(i) - centre) * array.spacing_m;
            // Conjugate steering weight exp(-j k x sin(steer)).
            const Real ph = -k * x * s;
            const Real cw = std::cos(ph);
            const Real sw = std::sin(ph);
            const Real er = elements[i].real();
            const Real ei = elements[i].imag();
            re += er * cw - ei * sw;
            im += er * sw + ei * cw;
        }
        out_power[a] = re * re + im * im;
    }
    return n_angles;
}

Real estimate_bearing_rad(std::span<const Real> power,
                          Real angle_min_rad, Real angle_max_rad) noexcept {
    if (power.empty() || !(angle_max_rad > angle_min_rad)) return kZero;
    const std::size_t n = power.size();
    if (n == 1) return angle_min_rad;

    std::size_t peak = 0;
    for (std::size_t i = 1; i < n; ++i) {
        if (power[i] > power[peak]) peak = i;
    }

    const Real step = (angle_max_rad - angle_min_rad) / static_cast<Real>(n - 1);
    Real offset = kZero;
    if (peak > 0 && peak + 1 < n) {
        // Parabola through the three magnitudes, as in the arrival-time
        // estimator: the mainlobe is far better fitted by a parabola in
        // amplitude than in power.
        const Real a = std::sqrt(power[peak - 1]);
        const Real b = std::sqrt(power[peak]);
        const Real c = std::sqrt(power[peak + 1]);
        const Real den = a - kTwo * b + c;
        if (den < kZero) {
            offset = clamp((a - c) / (kTwo * den), static_cast<Real>(-0.5),
                           static_cast<Real>(0.5));
        }
    }
    return angle_min_rad + step * (static_cast<Real>(peak) + offset);
}

Real split_beam_bearing_rad(const LineArray& array, Real lambda_m,
                            std::span<const Complex> elements,
                            Real steer_rad) noexcept {
    if (!array.valid() || !(lambda_m > kZero)) return steer_rad;
    if (elements.size() < array.element_count) return steer_rad;
    const std::size_t half = array.element_count / 2;
    if (half < 1) return steer_rad;

    const Real k = kTwo * kPi / lambda_m;
    const Real centre = static_cast<Real>(array.element_count - 1) / kTwo;
    const Real s0 = std::sin(steer_rad);

    Real lo_re = kZero, lo_im = kZero, hi_re = kZero, hi_im = kZero;
    for (std::size_t i = 0; i < array.element_count; ++i) {
        const Real x = (static_cast<Real>(i) - centre) * array.spacing_m;
        const Real ph = -k * x * s0;
        const Real cw = std::cos(ph);
        const Real sw = std::sin(ph);
        const Real er = elements[i].real();
        const Real ei = elements[i].imag();
        const Real re = er * cw - ei * sw;
        const Real im = er * sw + ei * cw;
        if (i < half) { lo_re += re; lo_im += im; } else { hi_re += re; hi_im += im; }
    }

    // arg(upper * conj(lower)) is the phase the wavefront accumulates between
    // the two half-array centroids.
    const Real cross_re = hi_re * lo_re + hi_im * lo_im;
    const Real cross_im = hi_im * lo_re - hi_re * lo_im;
    if (cross_re == kZero && cross_im == kZero) return steer_rad;
    const Real dphi = std::atan2(cross_im, cross_re);

    const Real separation = static_cast<Real>(array.element_count) * array.spacing_m / kTwo;
    const Real delta_sin = dphi / (k * separation);
    const Real target = s0 + delta_sin;
    if (!(std::fabs(target) <= kOne)) return steer_rad;
    return std::asin(target);
}

}  // namespace phantom
