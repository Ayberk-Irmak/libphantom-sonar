// SPDX-License-Identifier: Apache-2.0
// Line arrays, beamforming and bearing.
//
// The bearing bound is the spatial twin of the arrival-time bound verified in
// v0.2: the element index takes the place of time and Sum n'^2 = N(N^2-1)/12
// about the array centre takes the place of the waveform's mean-square
// bandwidth. Comparing an estimator against it is the same kind of check, and
// worth more than any amount of agreement with another implementation.
#include "framework.hpp"

#include "phantom/array.hpp"
#include "phantom/reverberation.hpp"

#include <algorithm>
#include <array>
#include <cmath>

using namespace phantom;

namespace {

constexpr double kPiD = 3.14159265358979323846;
constexpr Real   kC   = 1500;

std::array<Complex, 256> g_elements;
std::array<Real, 4096>   g_power;

LineArray half_wave(std::size_t n, Real lambda) {
    LineArray a;
    a.element_count = n;
    a.spacing_m = lambda / 2;
    return a;
}

}  // namespace

// ---------------------------------------------------------------------------
// Beam pattern
// ---------------------------------------------------------------------------

PT_TEST(array_factor_peaks_at_the_steer_angle) {
    const Real lambda = wavelength_m(10000, kC);
    const LineArray a = half_wave(32, lambda);
    PT_CHECK_REL(lambda, 0.15, pt::tol(1e-12, 1e-5));
    PT_CHECK_REL(a.aperture_m(), 31.0 * 0.075, pt::tol(1e-12, 1e-5));

    for (double steer_deg : {-60.0, -20.0, 0.0, 15.0, 45.0}) {
        const auto steer = deg2rad(static_cast<Real>(steer_deg));
        PT_CHECK_NEAR(array_factor(a, lambda, steer, steer), 1.0, pt::tol(1e-9, 1e-5));
        // And nothing anywhere else exceeds it.
        for (int look = -90; look <= 90; ++look) {
            PT_CHECK(array_factor(a, lambda, steer, deg2rad(static_cast<Real>(look)))
                     <= static_cast<Real>(1) + static_cast<Real>(1e-9));
        }
    }
}

PT_TEST(array_factor_nulls_are_where_the_closed_form_says) {
    // psi = 2 pi m / N, i.e. sin(look) - sin(steer) = m lambda / (N d).
    const Real lambda = wavelength_m(10000, kC);
    const LineArray a = half_wave(24, lambda);

    std::printf("       N = 24, d = lambda/2, broadside\n");
    std::printf("       %6s %14s %16s\n", "m", "null (deg)", "|B| there");
    for (int m = 1; m <= 4; ++m) {
        const double delta = static_cast<double>(m) * static_cast<double>(lambda)
                           / (24.0 * static_cast<double>(a.spacing_m));
        if (delta > 1.0) continue;
        const auto look = static_cast<Real>(std::asin(delta));
        const double b = static_cast<double>(array_factor(a, lambda, 0, look));
        std::printf("       %6d %14.4f %16.3e\n", m,
                    std::asin(delta) * 180.0 / kPiD, b);
        PT_CHECK(b < pt::tol(1e-9, 1e-4));
    }

    // The helper agrees with the m = 1 case, and does so for a steered beam too
    // -- where the offset is NOT lambda/(N d) but asin(sin(steer) + lambda/(N d))
    // minus the steer angle.
    for (double steer_deg : {0.0, 30.0, 60.0}) {
        const auto steer = deg2rad(static_cast<Real>(steer_deg));
        const Real off = first_null_offset_rad(a, lambda, steer);
        PT_CHECK(off > 0);
        const double b = static_cast<double>(array_factor(a, lambda, steer, steer + off));
        std::printf("       steer %4.0f deg -> first null %+.4f deg, |B| = %.3e\n",
                    steer_deg, static_cast<double>(rad2deg(off)), b);
        PT_CHECK(b < pt::tol(1e-9, 1e-4));
        // A steered beam is broader: the projected aperture shrinks as cos.
        if (steer_deg > 0) {
            PT_CHECK(off > first_null_offset_rad(a, lambda, 0));
        }
    }
}

PT_TEST(first_sidelobe_approaches_minus_thirteen_decibels) {
    // -13.26 dB is the signature of uniform shading. A window would lower it at
    // the cost of a wider mainlobe; getting a different number from a uniform
    // array means the pattern is wrong.
    const Real lambda = wavelength_m(10000, kC);
    std::printf("       %6s %18s\n", "N", "first sidelobe (dB)");
    double last = 0;
    for (std::size_t n : {std::size_t(8), std::size_t(16), std::size_t(32),
                          std::size_t(64), std::size_t(128)}) {
        const LineArray a = half_wave(n, lambda);
        // Search between the first and second nulls in psi.
        const double lo = 2.0 * kPiD / static_cast<double>(n);
        const double hi = 4.0 * kPiD / static_cast<double>(n);
        double best = 0;
        for (int i = 1; i < 20000; ++i) {
            const double psi = lo + (hi - lo) * static_cast<double>(i) / 20000.0;
            // psi = k d sin(look) at broadside, so sin(look) = psi / (k d).
            const double k = 2.0 * kPiD / static_cast<double>(lambda);
            const double s = psi / (k * static_cast<double>(a.spacing_m));
            if (std::fabs(s) > 1.0) continue;
            const auto look = static_cast<Real>(std::asin(s));
            best = std::max(best, static_cast<double>(array_factor(a, lambda, 0, look)));
        }
        last = 20.0 * std::log10(best);
        std::printf("       %6zu %18.3f\n", n, last);
        PT_CHECK(last < -12.5);
        PT_CHECK(last > -13.5);
    }
    // The limit is -13.26 dB and a 128-element array is essentially there.
    PT_CHECK_NEAR(last, -13.26, 0.05);
}

PT_TEST(grating_lobes_appear_past_half_wavelength_spacing) {
    const Real lambda = wavelength_m(10000, kC);

    // The classic rule: steering to endfire demands d <= lambda/2.
    PT_CHECK_REL(max_spacing_no_grating_lobes_m(lambda, kHalfPi),
                 static_cast<double>(lambda) / 2.0, pt::tol(1e-9, 1e-5));
    // A beam that never leaves broadside can use a full wavelength.
    PT_CHECK_REL(max_spacing_no_grating_lobes_m(lambda, 0),
                 static_cast<double>(lambda), pt::tol(1e-9, 1e-5));

    std::printf("       %10s %10s %14s\n", "d/lambda", "steer", "grating lobe");
    for (double ratio : {0.5, 0.9, 1.0, 1.5}) {
        LineArray a;
        a.element_count = 16;
        a.spacing_m = static_cast<Real>(ratio) * lambda;
        for (double steer_deg : {0.0, 60.0}) {
            const bool g = has_grating_lobe(a, lambda, deg2rad(static_cast<Real>(steer_deg)));
            std::printf("       %10.2f %9.0fd %14s\n", ratio, steer_deg, g ? "yes" : "no");
            // A grating lobe is a full-height repeat of the mainlobe, so when
            // the flag says yes there must actually be one.
            if (g) {
                double best_far = 0;
                for (int look = -90; look <= 90; ++look) {
                    const auto lk = deg2rad(static_cast<Real>(look));
                    if (std::fabs(static_cast<double>(look) - steer_deg) < 15.0) continue;
                    best_far = std::max(best_far,
                                        static_cast<double>(array_factor(a, lambda,
                                            deg2rad(static_cast<Real>(steer_deg)), lk)));
                }
                PT_CHECK(best_far > 0.9);
            }
        }
    }
    PT_CHECK(!has_grating_lobe(half_wave(16, lambda), lambda, 0));
}

// ---------------------------------------------------------------------------
// Gain
// ---------------------------------------------------------------------------

PT_TEST(array_gain_is_ten_log_n_and_is_actually_achieved) {
    // Signal adds coherently, noise does not. Verified by Monte Carlo rather
    // than asserted, because the formula is trivial and the implementation of
    // the beamformer is not.
    const Real lambda = wavelength_m(10000, kC);
    pt::Rng rng(20260807);

    std::printf("       %6s %12s %14s %12s\n", "N", "10log10(N)", "measured (dB)", "trials");
    for (std::size_t n : {std::size_t(4), std::size_t(16), std::size_t(64)}) {
        const LineArray a = half_wave(n, lambda);
        const double amplitude = 1.0;
        const double sigma = 4.0;                 // element SNR = -12 dB

        const std::size_t trials = 400;
        double sum_signal = 0, sum_noise = 0;
        for (std::size_t t = 0; t < trials; ++t) {
            synthesize_plane_wave(a, lambda, 0, static_cast<Real>(amplitude), g_elements);
            // Noise-only run first, using the same beamformer.
            for (std::size_t i = 0; i < n; ++i) {
                const auto re = static_cast<Real>(sigma * rng.normal() / std::sqrt(2.0));
                const auto im = static_cast<Real>(sigma * rng.normal() / std::sqrt(2.0));
                g_elements[i] = Complex(re, im);
            }
            beamform_power(a, lambda, g_elements, 0, static_cast<Real>(1e-6),
                           std::span<Real>(g_power.data(), 1));
            sum_noise += static_cast<double>(g_power[0]);

            synthesize_plane_wave(a, lambda, 0, static_cast<Real>(amplitude), g_elements);
            beamform_power(a, lambda, g_elements, 0, static_cast<Real>(1e-6),
                           std::span<Real>(g_power.data(), 1));
            sum_signal += static_cast<double>(g_power[0]);
        }
        const double out_snr = (sum_signal / static_cast<double>(trials))
                             / (sum_noise / static_cast<double>(trials));
        const double in_snr = amplitude * amplitude / (sigma * sigma);
        const double gain_db = 10.0 * std::log10(out_snr / in_snr);
        std::printf("       %6zu %12.3f %14.3f %12zu\n",
                    n, 10.0 * std::log10(static_cast<double>(n)), gain_db, trials);
        PT_CHECK_NEAR(gain_db, 10.0 * std::log10(static_cast<double>(n)), 0.6);
        PT_CHECK_NEAR(array_gain_db(a), 10.0 * std::log10(static_cast<double>(n)),
                      pt::tol(1e-9, 1e-4));
    }
}

// ---------------------------------------------------------------------------
// Bearing, against the bound
// ---------------------------------------------------------------------------

PT_TEST(bearing_crlb_scales_as_the_derivation_says) {
    // Two structural properties, both read straight off
    // var >= 6 / (rho (k d cos)^2 N(N^2-1)):
    //   accuracy improves as N^(-3/2), because elements buy signal AND aperture
    //   accuracy degrades as 1/cos(theta), because the projected aperture shrinks
    const Real lambda = wavelength_m(10000, kC);

    std::printf("       %6s %16s %14s\n", "N", "CRLB (deg)", "ratio to N/2");
    double prev = 0;
    for (std::size_t n : {std::size_t(8), std::size_t(16), std::size_t(32), std::size_t(64)}) {
        const double s = static_cast<double>(
            rad2deg(bearing_crlb_rad(half_wave(n, lambda), lambda, 0, 1)));
        std::printf("       %6zu %16.5f %14.4f\n", n, s, prev > 0 ? prev / s : 0.0);
        if (prev > 0) {
            // Doubling N should improve the bound by 2^1.5 = 2.828.
            PT_CHECK_NEAR(prev / s, std::pow(2.0, 1.5), 0.05);
        }
        prev = s;
    }

    std::printf("       %10s %16s %14s\n", "bearing", "CRLB (deg)", "1/cos scaling");
    const LineArray a = half_wave(32, lambda);
    const double at0 = static_cast<double>(rad2deg(bearing_crlb_rad(a, lambda, 0, 1)));
    for (double deg : {0.0, 30.0, 60.0, 75.0}) {
        const auto th = deg2rad(static_cast<Real>(deg));
        const double s = static_cast<double>(rad2deg(bearing_crlb_rad(a, lambda, th, 1)));
        std::printf("       %9.0fd %16.5f %14.4f\n", deg, s,
                    at0 / std::cos(deg * kPiD / 180.0));
        PT_CHECK_NEAR(s, at0 / std::cos(deg * kPiD / 180.0), pt::tol(1e-6, 1e-3));
    }

    // A 32-element half-wave array at 0 dB element SNR resolves bearing to
    // about a thirteenth of its own beamwidth.
    const double crlb = static_cast<double>(rad2deg(bearing_crlb_rad(a, lambda, 0, 1)));
    const double bw = static_cast<double>(rad2deg(beamwidth_3db_rad(a, lambda, 0)));
    std::printf("       beamwidth %.3f deg, CRLB %.4f deg -> 1/%.1f of a beamwidth\n",
                bw, crlb, bw / crlb);
    PT_CHECK(bw / crlb > 10.0);

    PT_CHECK(bearing_crlb_rad(a, lambda, kHalfPi, 1) == 0);   // endfire: no information
    PT_CHECK(bearing_crlb_rad(a, lambda, 0, 0) == 0);
}

PT_TEST(beamform_bearing_estimator_approaches_the_bound) {
    // The same comparison as the arrival-time estimator in v0.2, in the spatial
    // domain: measure the estimator's scatter and hold it against the best any
    // unbiased estimator could do.
    const Real lambda = wavelength_m(10000, kC);
    const LineArray a = half_wave(32, lambda);
    const Real truth = deg2rad(static_cast<Real>(7));

    // Scan the mainlobe finely enough that the grid is not the limit.
    const Real scan_lo = truth - deg2rad(static_cast<Real>(6));
    const Real scan_hi = truth + deg2rad(static_cast<Real>(6));
    constexpr std::size_t kAngles = 2001;

    pt::Rng rng(777);
    std::printf("       %10s %14s %14s %8s\n", "sigma", "measured deg", "CRLB deg", "ratio");
    for (double sigma : {2.0, 1.0, 0.4}) {
        const double amplitude = 1.0;
        const double element_snr = amplitude * amplitude / (sigma * sigma);
        const double crlb = static_cast<double>(
            rad2deg(bearing_crlb_rad(a, lambda, truth, static_cast<Real>(element_snr))));

        const std::size_t trials = 400;
        double sum = 0, sum_sq = 0;
        std::size_t used = 0;
        for (std::size_t t = 0; t < trials; ++t) {
            synthesize_plane_wave(a, lambda, truth, static_cast<Real>(amplitude), g_elements);
            for (std::size_t i = 0; i < a.element_count; ++i) {
                const auto re = static_cast<Real>(sigma * rng.normal() / std::sqrt(2.0));
                const auto im = static_cast<Real>(sigma * rng.normal() / std::sqrt(2.0));
                g_elements[i] += Complex(re, im);
            }
            beamform_power(a, lambda, g_elements, scan_lo, scan_hi,
                           std::span<Real>(g_power.data(), kAngles));
            const double est = static_cast<double>(rad2deg(estimate_bearing_rad(
                std::span<const Real>(g_power.data(), kAngles), scan_lo, scan_hi)));
            const double err = est - 7.0;
            // Outside the mainlobe the estimate is a sidelobe capture, a
            // different failure mode from estimator scatter; counted, not hidden.
            if (std::fabs(err) > 3.0) continue;
            sum += err;
            sum_sq += err * err;
            ++used;
        }
        PT_CHECK(used > trials * 3 / 4);
        const double mean = sum / static_cast<double>(used);
        const double sd = std::sqrt(sum_sq / static_cast<double>(used) - mean * mean);
        std::printf("       %10.2f %14.5f %14.5f %8.2f   (%zu/%zu, bias %+.4f)\n",
                    sigma, sd, crlb, sd / crlb, used, trials, mean);
        PT_CHECK(std::fabs(mean) < 1.5 * sd);
        PT_CHECK(sd > 0.75 * crlb);      // cannot beat its own bound
        PT_CHECK(sd < 1.6 * crlb);       // and is close to it
    }
}

PT_TEST(split_beam_matches_the_scan_and_is_far_cheaper) {
    // Split-beam reads bearing off the phase between the two half-arrays. It is
    // unambiguous exactly out to the first null, because |dphi| = pi there.
    const Real lambda = wavelength_m(10000, kC);
    const LineArray a = half_wave(32, lambda);
    const Real steer = 0;

    std::printf("       %12s %14s %14s\n", "truth (deg)", "split-beam", "error (deg)");
    const double null_deg = static_cast<double>(rad2deg(first_null_offset_rad(a, lambda, steer)));
    for (double deg : {-3.0, -1.0, 0.0, 0.5, 2.0, 3.5}) {
        if (std::fabs(deg) > null_deg) continue;
        const auto truth = deg2rad(static_cast<Real>(deg));
        synthesize_plane_wave(a, lambda, truth, 1, g_elements);
        const double est = static_cast<double>(
            rad2deg(split_beam_bearing_rad(a, lambda, g_elements, steer)));
        std::printf("       %12.2f %14.5f %14.2e\n", deg, est, std::fabs(est - deg));
        PT_CHECK_NEAR(est, deg, pt::tol(1e-6, 1e-3));
    }
    std::printf("       unambiguous out to the first null at %.3f deg\n", null_deg);

    // Beyond the first null the phase wraps and the answer is wrong -- which is
    // the documented limit, not a defect. A target there must be found by the
    // scan first.
    const auto far = deg2rad(static_cast<Real>(null_deg * 2.5));
    synthesize_plane_wave(a, lambda, far, 1, g_elements);
    const double wrapped = static_cast<double>(
        rad2deg(split_beam_bearing_rad(a, lambda, g_elements, steer)));
    std::printf("       a source at %.2f deg reads as %.2f deg (wrapped)\n",
                null_deg * 2.5, wrapped);
    PT_CHECK(std::fabs(wrapped - null_deg * 2.5) > 0.5);

    // Rejections.
    PT_CHECK(split_beam_bearing_rad(a, lambda, std::span<const Complex>(), 0) == 0);
    LineArray bad;
    PT_CHECK(split_beam_bearing_rad(bad, lambda, g_elements, 0) == 0);
}

// ---------------------------------------------------------------------------
// Why the array was worth building
// ---------------------------------------------------------------------------

PT_TEST(narrow_beams_buy_exactly_what_reverberation_charges) {
    // v0.6 showed reverberation scales with the ensonified area, and that a
    // tenth of the beamwidth is worth 10 dB. An array is how a beamwidth gets
    // narrower, so the two results have to join up.
    const Real lambda = wavelength_m(10000, kC);
    const Real ts = 10, ss = -30, tau = static_cast<Real>(0.001), range = 2000;

    std::printf("       %6s %14s %16s %14s\n",
                "N", "beamwidth", "E/R ratio (dB)", "array gain");
    double prev_er = 0;
    for (std::size_t n : {std::size_t(8), std::size_t(80)}) {
        const LineArray a = half_wave(n, lambda);
        const Real bw = beamwidth_3db_rad(a, lambda, 0);
        const Real area = ensonified_area_m2(range, bw, tau, kC);
        const double er = static_cast<double>(echo_to_reverberation_ratio_db(ts, ss, area));
        std::printf("       %6zu %13.3fd %16.2f %13.2f\n",
                    n, static_cast<double>(rad2deg(bw)), er,
                    static_cast<double>(array_gain_db(a)));
        if (prev_er != 0) {
            // Ten times the elements is ten times the aperture is a tenth of
            // the beamwidth is exactly 10 dB of echo-to-reverberation ratio.
            PT_CHECK_NEAR(er - prev_er, 10.0, 0.05);
        }
        prev_er = er;
    }

    std::printf("       A ten-fold array buys 10 dB against reverberation from the\n"
                "       beamwidth and %.1f dB against isotropic noise from the gain --\n"
                "       two different mechanisms that happen to agree here.\n",
                static_cast<double>(array_gain_db(half_wave(80, lambda)))
              - static_cast<double>(array_gain_db(half_wave(8, lambda))));
}

PT_TEST(phase_steering_has_a_bandwidth_limit_worth_knowing) {
    // The library's own waveforms are 12 kHz chirps. A phase-steered array
    // cannot handle them at a steered angle, and saying so is cheaper than
    // shipping a beamformer that quietly smears the beam.
    const Real lambda = wavelength_m(10000, kC);
    const LineArray a = half_wave(32, lambda);

    std::printf("       %10s %16s %18s\n", "steer", "traversal (ms)", "usable BW (Hz)");
    for (double deg : {0.0, 15.0, 45.0, 90.0}) {
        const auto th = deg2rad(static_cast<Real>(deg));
        const double trav = static_cast<double>(a.aperture_m())
                          * std::fabs(std::sin(deg * kPiD / 180.0)) / 1500.0;
        const double bw = static_cast<double>(narrowband_bandwidth_limit_hz(a, th, kC));
        std::printf("       %9.0fd %16.4f %18.1f\n", deg, trav * 1e3, bw);
        if (deg > 0) {
            PT_CHECK(bw > 0);
            PT_CHECK_NEAR(bw, 1.0 / (10.0 * trav), pt::tol(1e-6, 1e-3));
        }
    }
    // Broadside costs nothing: there is no traversal to decorrelate over.
    PT_CHECK(narrowband_bandwidth_limit_hz(a, 0, kC) == 0);
    // At 45 degrees the limit is far below the 12 kHz the analyser transmits.
    PT_CHECK(narrowband_bandwidth_limit_hz(a, deg2rad(static_cast<Real>(45)), kC) < 12000);
}
