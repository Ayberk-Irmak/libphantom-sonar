// SPDX-License-Identifier: Apache-2.0
// Wideband beamforming, shading, and adaptive beams.
//
// The wideband test is the one that matters: v0.7 measured that phase steering
// dies at 91 Hz of bandwidth off broadside, and this library transmits 12 kHz
// chirps. If per-bin steering does not fix that, the array cannot be used with
// the library's own waveforms and the whole module is decorative.
#include "framework.hpp"

#include "phantom/beamformer.hpp"
#include "phantom/matched_filter.hpp"
#include "phantom/ping_analyzer.hpp"

#include <array>
#include <cmath>

using namespace phantom;

namespace {

constexpr double      kPiD = 3.14159265358979323846;
constexpr Real        kC   = 1500;
constexpr Real        kFs  = 96000;
constexpr std::size_t kFft = 2048;
constexpr std::size_t kElem = 16;

FftPlan<kFft>                          g_plan;
std::array<Real, kElem * kFft>         g_elements;
std::array<Complex, kElem * kFft>      g_spectra;
std::array<Complex, kFft>              g_work;
std::array<Real, kFft>                 g_beam;
std::array<Complex, 256>               g_snapshot;
std::array<Complex, kElem * kElem>     g_cov;
std::array<Complex, kElem * kElem + 2 * kElem> g_mvdr_work;
std::array<Real, 4096>                 g_power;

LineArray half_wave(std::size_t n, Real lambda) {
    LineArray a;
    a.element_count = n;
    a.spacing_m = lambda / 2;
    return a;
}

// Lays a plane wave from `bearing` across the elements, evaluating the
// waveform's closed-form phase at each element's own arrival time -- so the
// fractional delays are exact rather than interpolated.
//
// SIGN CONVENTION. The element at +x sees the wavefront EARLIER, so it is
// further into the waveform at a given sample index: t = (i - at)/fs + x
// sin(theta)/c. That matches synthesize_plane_wave's exp(+j k x sin(theta)),
// which v0.7 already validated against the Cramer-Rao bound and split-beam.
// Getting it backwards steers every beam to the mirror bearing, which looks
// like a working beamformer until you check where it points.
void illuminate(const LineArray& array, const PulseSpec& spec, Real bearing,
                std::size_t at, Real amplitude) {
    const Real centre = static_cast<Real>(array.element_count - 1) / 2;
    for (Real& v : g_elements) v = 0;
    for (std::size_t n = 0; n < array.element_count; ++n) {
        const Real x = (static_cast<Real>(n) - centre) * array.spacing_m;
        const Real advance = x * std::sin(bearing) / kC * kFs;   // samples, early
        for (std::size_t i = 0; i < kFft; ++i) {
            const Real t = (static_cast<Real>(i) - static_cast<Real>(at) + advance) / kFs;
            if (t < 0 || t > spec.duration_s) continue;
            g_elements[n * kFft + i] = amplitude * std::cos(pulse_phase(spec, t));
        }
    }
}

PulseSpec chirp(Real f0, Real f1, Real dur) {
    PulseSpec s;
    s.type = (f1 > f0) ? PulseType::LfmUp : PulseType::LfmDown;
    s.f_start_hz = f0;
    s.f_end_hz = f1;
    s.duration_s = dur;
    return s;
}

}  // namespace

// ---------------------------------------------------------------------------
// Shading
// ---------------------------------------------------------------------------

PT_TEST(shading_buys_sidelobes_and_pays_in_beamwidth) {
    // The shaded pattern is the DTFT of the window, so its sidelobe level IS
    // the window's -- these are the familiar numbers, not array-specific ones.
    const Real lambda = wavelength_m(10000, kC);
    const LineArray a = half_wave(64, lambda);

    struct Expect { Shading s; double sidelobe_db; double tol; };
    const Expect cases[] = {
        {Shading::Uniform,  -13.26, 0.3},
        {Shading::Hann,     -31.5,  1.5},
        {Shading::Hamming,  -42.7,  3.0},
        {Shading::Blackman, -58.1,  4.0},
    };

    std::printf("       %10s %14s %12s %14s %12s\n",
                "window", "sidelobe (dB)", "published", "width x uniform", "loss (dB)");
    double uniform_width = 0;
    for (const Expect& e : cases) {
        // Peak sidelobe: scan outside the mainlobe, which is bounded by the
        // first null of the widest window here.
        double peak_sl = 0;
        double first_null = 0;
        double prev = 1.0;
        for (int i = 1; i <= 40000; ++i) {
            const double deg = 0.001 * static_cast<double>(i);
            if (deg > 20.0) break;
            const auto look = deg2rad(static_cast<Real>(deg));
            const double b = static_cast<double>(shaded_array_factor(a, e.s, lambda, 0, look));
            if (first_null == 0 && b > prev && prev < 0.05) first_null = deg;
            if (first_null > 0) peak_sl = std::max(peak_sl, b);
            prev = b;
        }
        // -3 dB half-width, for the beamwidth ratio.
        double half = 0;
        for (int i = 1; i <= 40000; ++i) {
            const double deg = 0.001 * static_cast<double>(i);
            const auto look = deg2rad(static_cast<Real>(deg));
            if (static_cast<double>(shaded_array_factor(a, e.s, lambda, 0, look))
                < 1.0 / std::sqrt(2.0)) { half = deg; break; }
        }
        if (e.s == Shading::Uniform) uniform_width = half;

        const double sl_db = 20.0 * std::log10(peak_sl);
        std::printf("       %10s %14.2f %12.1f %14.3f %12.2f\n",
                    shading_name(e.s), sl_db, e.sidelobe_db, half / uniform_width,
                    static_cast<double>(shading_loss_db(e.s, 64)));
        PT_CHECK_NEAR(sl_db, e.sidelobe_db, e.tol);

        // Every window is a trade: wider mainlobe, and gain given up.
        if (e.s != Shading::Uniform) {
            PT_CHECK(half > uniform_width);
            PT_CHECK(static_cast<double>(shading_loss_db(e.s, 64)) < 0);
        } else {
            PT_CHECK_NEAR(shading_loss_db(e.s, 64), 0.0, 1e-9);
        }
    }

    PT_CHECK_NEAR(shading_weight(Shading::Uniform, 3, 16), 1.0, 1e-12);
    PT_CHECK_NEAR(shading_weight(Shading::Hann, 0, 16), 0.0, 1e-12);
    PT_CHECK(static_cast<double>(shading_weight(Shading::Hann, 8, 17)) > 0.99);
    PT_CHECK(shading_weight(Shading::Hann, 99, 16) == 0);
}

// ---------------------------------------------------------------------------
// Wideband beamforming -- the point of the release
// ---------------------------------------------------------------------------

PT_TEST(wideband_path_reproduces_the_narrowband_beam_pattern) {
    // Ties the two halves of the library together: a narrowband tone pushed
    // through the wideband beamformer must trace out the same pattern that
    // array_factor gives in closed form. A sign error in either steering
    // convention shows up here as a mirrored pattern rather than as a beam
    // that merely points somewhere plausible.
    const Real f0 = 12000;
    const Real lambda = wavelength_m(f0, kC);
    const LineArray a = half_wave(kElem, lambda);
    const Real bearing = deg2rad(static_cast<Real>(18));

    PulseSpec tone;
    tone.type = PulseType::Cw;
    tone.f_start_hz = f0;
    tone.duration_s = static_cast<Real>(kFft) / kFs;   // fills the block
    illuminate(a, tone, bearing, 0, 1);
    PT_CHECK(prepare_element_spectra(g_plan.view(), a.element_count,
                                     g_elements, g_spectra, g_work));

    std::printf("       %10s %16s %16s %10s\n",
                "steer", "wideband |B|", "array_factor", "diff");
    double worst = 0;
    for (double deg = -40; deg <= 40; deg += 4) {
        beamform_wideband(g_plan.view(), a, g_spectra, kFs, kC,
                          deg2rad(static_cast<Real>(deg)), Shading::Uniform,
                          g_work, g_beam);
        // RMS of the beam, over the interior to avoid the block edges where
        // the tone was truncated.
        double e = 0;
        for (std::size_t i = kFft / 4; i < 3 * kFft / 4; ++i) {
            e += static_cast<double>(g_beam[i]) * static_cast<double>(g_beam[i]);
        }
        const double rms = std::sqrt(e / static_cast<double>(kFft / 2));
        // The tone has unit amplitude, so a beam on target has RMS 1/sqrt(2).
        const double measured = rms * std::sqrt(2.0);
        const double expect = static_cast<double>(
            array_factor(a, lambda, deg2rad(static_cast<Real>(deg)), bearing));
        if (std::fabs(deg - 18.0) < 13.0 || std::fabs(deg) < 1) {
            std::printf("       %9.0fd %16.5f %16.5f %10.2e\n",
                        deg, measured, expect, std::fabs(measured - expect));
        }
        worst = std::max(worst, std::fabs(measured - expect));
    }
    std::printf("       worst departure from the closed-form pattern: %.3e\n", worst);
    PT_CHECK(worst < pt::tol(0.02, 0.05));
}

PT_TEST(wideband_steering_works_where_phase_steering_fails) {
    // v0.7 measured the phase-steering bandwidth limit for this array at 45
    // degrees. Here is a signal far outside it, steered correctly.
    const Real lambda = wavelength_m(14000, kC);   // centre of the chirp band
    const LineArray a = half_wave(kElem, lambda);
    const PulseSpec s = chirp(8000, 20000, static_cast<Real>(0.005));
    const Real bearing = deg2rad(static_cast<Real>(35));

    const double limit = static_cast<double>(
        narrowband_bandwidth_limit_hz(a, bearing, kC));
    std::printf("       array phase-steering limit at 35 deg: %.0f Hz\n", limit);
    std::printf("       signal bandwidth: %.0f Hz -- %.0fx over\n",
                static_cast<double>(s.bandwidth_hz()),
                static_cast<double>(s.bandwidth_hz()) / limit);
    PT_CHECK(static_cast<double>(s.bandwidth_hz()) > 10.0 * limit);

    illuminate(a, s, bearing, 400, 1);
    PT_CHECK(prepare_element_spectra(g_plan.view(), a.element_count,
                                     g_elements, g_spectra, g_work));

    // Beam power on target versus off, as a function of steer angle.
    std::printf("       %10s %16s\n", "steer", "beam energy (dB)");
    double on_target = 0;
    double worst_off = -1e9;
    for (double deg : {0.0, 15.0, 30.0, 35.0, 40.0, 55.0}) {
        PT_CHECK(beamform_wideband(g_plan.view(), a, g_spectra, kFs, kC,
                                   deg2rad(static_cast<Real>(deg)),
                                   Shading::Uniform, g_work, g_beam));
        double energy = 0;
        for (std::size_t i = 0; i < kFft; ++i) {
            energy += static_cast<double>(g_beam[i]) * static_cast<double>(g_beam[i]);
        }
        const double db = 10.0 * std::log10(energy);
        std::printf("       %9.0fd %16.2f\n", deg, db);
        if (deg == 35.0) on_target = db;
        else if (std::fabs(deg - 35.0) > 10.0) worst_off = std::max(worst_off, db);
    }
    std::printf("       on target %.2f dB, best off-target %.2f dB -> %.1f dB contrast\n",
                on_target, worst_off, on_target - worst_off);
    PT_CHECK(on_target > worst_off + 6.0);

    // The beam steered on target must reproduce the waveform, not a smeared
    // version of it: correlate it against the replica and check the peak.
    static FftPlan<kFft> plan2;
    static std::array<Complex, kFft> rep{}, spec{}, corr{};
    beamform_wideband(g_plan.view(), a, g_spectra, kFs, kC, bearing,
                      Shading::Uniform, g_work, g_beam);
    const std::size_t l = render_analytic(s, kFs, rep);
    MatchedFilter mf;
    matched_filter_prepare(plan2.view(), std::span<const Complex>(rep.data(), l), spec, mf);
    const std::size_t lags = matched_filter_apply(plan2.view(), mf, g_beam, g_work, corr);
    std::size_t pk = 0;
    double pv = 0;
    for (std::size_t i = 0; i < lags; ++i) {
        const double v = static_cast<double>(std::abs(corr[i]));
        if (v > pv) { pv = v; pk = i; }
    }
    std::printf("       matched filter peak at lag %zu (pulse placed at 400)\n", pk);
    // A correctly steered wideband beam keeps the pulse compressible; a
    // phase-steered one would smear it and move the peak.
    PT_CHECK(pk + 3 >= 400 && pk <= 403);
}

PT_TEST(wideband_beam_carries_a_bearing_into_the_analyser) {
    // The integration this release exists for: run the pulse analyser on each
    // beam, and a detection now has a direction.
    const Real lambda = wavelength_m(14000, kC);
    const LineArray a = half_wave(kElem, lambda);
    const PulseSpec s = chirp(8000, 20000, static_cast<Real>(0.005));
    const Real truth = deg2rad(static_cast<Real>(-22));

    illuminate(a, s, truth, 500, 1);
    pt::Rng rng(4242);
    for (Real& v : g_elements) v += static_cast<Real>(0.5 * rng.normal());
    prepare_element_spectra(g_plan.view(), a.element_count, g_elements, g_spectra, g_work);

    static PulseBank<4, kFft> bank(kFs);
    bank.clear();
    bank.add(s);
    static AnalyzerScratch<kFft> scratch;
    std::array<PulseDescriptor, 8> pdw{};

    DetectorConfig cfg;
    cfg.cfar_guard = suggested_cfar_guard(bank.view());
    cfg.cfar_train = 128;
    cfg.threshold_alpha = cfar_alpha(256, static_cast<Real>(1e-4));
    cfg.dead_time_s = suggested_dead_time_s(bank.view());

    std::printf("       %10s %10s %12s %12s\n", "beam", "detections", "ToA (ms)", "SNR dB");
    double best_snr = -1e9;
    double best_beam = 0;
    for (double deg = -40; deg <= 40; deg += 5) {
        beamform_wideband(g_plan.view(), a, g_spectra, kFs, kC,
                          deg2rad(static_cast<Real>(deg)), Shading::Uniform,
                          g_work, g_beam);
        const std::size_t n = analyze_block(bank.view(), cfg, g_beam, 0,
                                            scratch.view(), pdw);
        if (n > 0 && static_cast<double>(pdw[0].snr_db) > best_snr) {
            best_snr = static_cast<double>(pdw[0].snr_db);
            best_beam = deg;
        }
        if (std::fabs(deg + 22.0) < 12.0 || std::fabs(deg) < 1.0) {
            std::printf("       %9.0fd %10zu %12.4f %12.2f\n", deg, n,
                        n > 0 ? static_cast<double>(pdw[0].toa_s) * 1e3 : 0.0,
                        n > 0 ? static_cast<double>(pdw[0].snr_db) : 0.0);
        }
    }
    std::printf("       strongest beam %+.0f deg (truth %.0f), SNR %.2f dB\n",
                best_beam, -22.0, best_snr);
    // Within one beam of the truth: the scan is 5 degrees and the beamwidth is
    // about 6.6, so the strongest beam is the adjacent one at worst.
    PT_CHECK(std::fabs(best_beam + 22.0) <= 5.0);
    PT_CHECK(best_snr > 10.0);
}

// ---------------------------------------------------------------------------
// Adaptive beamforming
// ---------------------------------------------------------------------------

PT_TEST(mvdr_agrees_with_conventional_on_a_single_source) {
    // Sanity before the interesting claim: with one source and plenty of
    // snapshots, both beamformers must put it in the same place.
    const Real lambda = wavelength_m(10000, kC);
    const LineArray a = half_wave(kElem, lambda);
    const Real truth = deg2rad(static_cast<Real>(11));

    pt::Rng rng(99);
    covariance_clear(kElem, g_cov);
    const std::size_t snaps = 200;
    for (std::size_t t = 0; t < snaps; ++t) {
        synthesize_plane_wave(a, lambda, truth, 1, g_snapshot);
        // A random source phase per snapshot, or the covariance is rank-1 and
        // degenerate in a way a real one never is.
        const Real ph = static_cast<Real>(2.0 * kPiD * rng.uniform01());
        const Complex rot(std::cos(ph), std::sin(ph));
        for (std::size_t i = 0; i < kElem; ++i) {
            g_snapshot[i] = g_snapshot[i] * rot
                          + Complex(static_cast<Real>(0.2 * rng.normal()),
                                    static_cast<Real>(0.2 * rng.normal()));
        }
        covariance_accumulate(std::span<const Complex>(g_snapshot.data(), kElem),
                              kElem, g_cov);
    }
    covariance_normalise(kElem, snaps, g_cov);

    constexpr std::size_t kAngles = 1801;
    const Real lo = deg2rad(static_cast<Real>(-45));
    const Real hi = deg2rad(static_cast<Real>(45));
    PT_CHECK(mvdr_power(a, lambda, g_cov, static_cast<Real>(0.01), lo, hi,
                        g_mvdr_work, std::span<Real>(g_power.data(), kAngles)) == kAngles);
    const double mvdr_deg = static_cast<double>(rad2deg(estimate_bearing_rad(
        std::span<const Real>(g_power.data(), kAngles), lo, hi)));

    synthesize_plane_wave(a, lambda, truth, 1, g_snapshot);
    beamform_power(a, lambda, g_snapshot, lo, hi,
                   std::span<Real>(g_power.data(), kAngles));
    const double conv_deg = static_cast<double>(rad2deg(estimate_bearing_rad(
        std::span<const Real>(g_power.data(), kAngles), lo, hi)));

    std::printf("       truth 11.00 deg: conventional %.3f, MVDR %.3f\n",
                conv_deg, mvdr_deg);
    PT_CHECK_NEAR(conv_deg, 11.0, 0.2);
    PT_CHECK_NEAR(mvdr_deg, 11.0, 0.3);
}

PT_TEST(mvdr_resolves_what_the_aperture_cannot) {
    // The claim worth making. Conventional beamforming cannot separate two
    // sources closer than about one null-to-peak spacing, no matter the SNR --
    // the resolution is set by the aperture and nothing else. MVDR is not bound
    // by that, because it places nulls rather than scanning a fixed pattern.
    const Real lambda = wavelength_m(10000, kC);
    const LineArray a = half_wave(kElem, lambda);

    const double limit_deg = static_cast<double>(
        rad2deg(conventional_resolution_limit_rad(a, lambda, 0)));
    const double sep_deg = limit_deg * 0.6;      // deliberately inside the limit
    const Real b1 = deg2rad(static_cast<Real>(-sep_deg / 2));
    const Real b2 = deg2rad(static_cast<Real>(+sep_deg / 2));
    std::printf("       conventional resolution limit %.2f deg; sources %.2f deg apart\n",
                limit_deg, sep_deg);

    pt::Rng rng(1234);
    covariance_clear(kElem, g_cov);
    const std::size_t snaps = 400;
    static std::array<Complex, 256> tmp{};
    for (std::size_t t = 0; t < snaps; ++t) {
        synthesize_plane_wave(a, lambda, b1, 1, g_snapshot);
        synthesize_plane_wave(a, lambda, b2, 1, tmp);
        // Independent phases: two genuinely uncorrelated sources. Coherent ones
        // defeat MVDR too, and would need spatial smoothing.
        const Real p1 = static_cast<Real>(2.0 * kPiD * rng.uniform01());
        const Real p2 = static_cast<Real>(2.0 * kPiD * rng.uniform01());
        const Complex r1(std::cos(p1), std::sin(p1));
        const Complex r2(std::cos(p2), std::sin(p2));
        for (std::size_t i = 0; i < kElem; ++i) {
            g_snapshot[i] = g_snapshot[i] * r1 + tmp[i] * r2
                          + Complex(static_cast<Real>(0.05 * rng.normal()),
                                    static_cast<Real>(0.05 * rng.normal()));
        }
        covariance_accumulate(std::span<const Complex>(g_snapshot.data(), kElem),
                              kElem, g_cov);
    }
    covariance_normalise(kElem, snaps, g_cov);

    constexpr std::size_t kAngles = 2001;
    const Real lo = deg2rad(static_cast<Real>(-15));
    const Real hi = deg2rad(static_cast<Real>(15));

    auto count_peaks = [&](std::span<const Real> p) {
        std::size_t peaks = 0;
        double top = 0;
        for (const Real v : p) top = std::max(top, static_cast<double>(v));
        for (std::size_t i = 1; i + 1 < p.size(); ++i) {
            if (p[i] > p[i - 1] && p[i] > p[i + 1]
                && static_cast<double>(p[i]) > top * 0.25) {
                ++peaks;
            }
        }
        return peaks;
    };

    // Conventional, using the covariance so the comparison is like for like:
    // its output is a^H R a.
    for (std::size_t ai = 0; ai < kAngles; ++ai) {
        const Real steer = lo + (hi - lo) * static_cast<Real>(ai)
                              / static_cast<Real>(kAngles - 1);
        synthesize_plane_wave(a, lambda, steer, 1, tmp);
        Complex acc(0, 0);
        Real total = 0;
        for (std::size_t i = 0; i < kElem; ++i) {
            acc = Complex(0, 0);
            for (std::size_t j = 0; j < kElem; ++j) acc += g_cov[i * kElem + j] * tmp[j];
            total += (std::conj(tmp[i]) * acc).real();
        }
        g_power[ai] = total;
    }
    const std::size_t conv_peaks = count_peaks(std::span<const Real>(g_power.data(), kAngles));

    static std::array<Real, 4096> mvdr_p{};
    mvdr_power(a, lambda, g_cov, static_cast<Real>(0.001), lo, hi, g_mvdr_work,
               std::span<Real>(mvdr_p.data(), kAngles));
    const std::size_t mvdr_peaks = count_peaks(std::span<const Real>(mvdr_p.data(), kAngles));

    std::printf("       conventional finds %zu peak(s), MVDR finds %zu\n",
                conv_peaks, mvdr_peaks);
    PT_CHECK(conv_peaks == 1);    // the aperture limit, exactly as advertised
    PT_CHECK(mvdr_peaks == 2);    // and MVDR is not bound by it

    // MVDR's two peaks must land on the two sources.
    double found[2] = {0, 0};
    std::size_t k = 0;
    double top = 0;
    for (std::size_t i = 0; i < kAngles; ++i) top = std::max(top, static_cast<double>(mvdr_p[i]));
    for (std::size_t i = 1; i + 1 < kAngles && k < 2; ++i) {
        if (mvdr_p[i] > mvdr_p[i - 1] && mvdr_p[i] > mvdr_p[i + 1]
            && static_cast<double>(mvdr_p[i]) > top * 0.25) {
            const double deg = -15.0 + 30.0 * static_cast<double>(i)
                                     / static_cast<double>(kAngles - 1);
            found[k++] = deg;
        }
    }
    std::printf("       MVDR peaks at %+.3f and %+.3f deg (truth %+.3f, %+.3f)\n",
                found[0], found[1], -sep_deg / 2, sep_deg / 2);
    PT_CHECK_NEAR(found[0], -sep_deg / 2, 0.25);
    PT_CHECK_NEAR(found[1], sep_deg / 2, 0.25);
}

PT_TEST(mvdr_needs_its_diagonal_loading) {
    // Loading is not a refinement. With fewer snapshots than elements the
    // covariance is singular and the Cholesky fails outright; the function
    // reports that rather than returning noise shaped like a spectrum.
    const Real lambda = wavelength_m(10000, kC);
    const LineArray a = half_wave(kElem, lambda);

    covariance_clear(kElem, g_cov);
    synthesize_plane_wave(a, lambda, 0, 1, g_snapshot);
    covariance_accumulate(std::span<const Complex>(g_snapshot.data(), kElem), kElem, g_cov);
    covariance_normalise(kElem, 1, g_cov);   // rank 1, 16 elements

    const Real lo = deg2rad(static_cast<Real>(-30));
    const Real hi = deg2rad(static_cast<Real>(30));
    PT_CHECK(mvdr_power(a, lambda, g_cov, 0, lo, hi, g_mvdr_work,
                        std::span<Real>(g_power.data(), 101)) == 0);
    PT_CHECK(mvdr_power(a, lambda, g_cov, static_cast<Real>(0.01), lo, hi, g_mvdr_work,
                        std::span<Real>(g_power.data(), 101)) == 101);
    std::printf("       rank-1 covariance over %zu elements: unloaded fails, loaded works\n",
                kElem);

    // Rejections.
    PT_CHECK(mvdr_power(a, lambda, g_cov, static_cast<Real>(0.01), lo, hi,
                        std::span<Complex>(g_mvdr_work.data(), 4),
                        std::span<Real>(g_power.data(), 101)) == 0);
    PT_CHECK(mvdr_power(a, lambda, g_cov, static_cast<Real>(0.01), hi, lo, g_mvdr_work,
                        std::span<Real>(g_power.data(), 101)) == 0);
    PT_CHECK(!covariance_accumulate(std::span<const Complex>(g_snapshot.data(), 2),
                                    kElem, g_cov));
}

PT_TEST(beamformer_rejects_bad_buffers) {
    const Real lambda = wavelength_m(10000, kC);
    const LineArray a = half_wave(kElem, lambda);

    PT_CHECK(!prepare_element_spectra(g_plan.view(), kElem,
                                      std::span<const Real>(g_elements.data(), 16),
                                      g_spectra, g_work));
    PT_CHECK(!prepare_element_spectra(g_plan.view(), kElem, g_elements,
                                      std::span<Complex>(g_spectra.data(), 16), g_work));
    PT_CHECK(!beamform_wideband(g_plan.view(), a, g_spectra, 0, kC, 0,
                                Shading::Uniform, g_work, g_beam));
    PT_CHECK(!beamform_wideband(g_plan.view(), a, g_spectra, kFs, 0, 0,
                                Shading::Uniform, g_work, g_beam));
    PT_CHECK(!beamform_wideband(g_plan.view(), a,
                                std::span<const Complex>(g_spectra.data(), 16),
                                kFs, kC, 0, Shading::Uniform, g_work, g_beam));
    LineArray bad;
    PT_CHECK(!beamform_wideband(g_plan.view(), bad, g_spectra, kFs, kC, 0,
                                Shading::Uniform, g_work, g_beam));

    PT_CHECK(conventional_resolution_limit_rad(bad, lambda, 0) == 0);
    PT_CHECK(static_cast<double>(conventional_resolution_limit_rad(a, lambda, kHalfPi)) > 1.0);
}

// ---------------------------------------------------------------------------
// v0.10: spatial smoothing, so MVDR survives coherent multipath
// ---------------------------------------------------------------------------

PT_TEST(mvdr_fails_on_coherent_sources_and_smoothing_repairs_it) {
    // MVDR's weakness, and the reason it is not simply better than conventional
    // beamforming. Two COHERENT arrivals -- one signal reaching the array by two
    // paths, which is exactly the multipath v0.4 produces -- make the covariance
    // rank-deficient in a way diagonal loading cannot fix. The adaptive
    // beamformer then nulls the target along with its own echo.
    const Real lambda = wavelength_m(10000, kC);
    const LineArray a = half_wave(kElem, lambda);
    const Real b1 = deg2rad(static_cast<Real>(-6));
    const Real b2 = deg2rad(static_cast<Real>(+6));

    pt::Rng rng(60606);
    static std::array<Complex, 256> tmp{};
    covariance_clear(kElem, g_cov);
    const std::size_t snaps = 400;
    for (std::size_t t = 0; t < snaps; ++t) {
        synthesize_plane_wave(a, lambda, b1, 1, g_snapshot);
        synthesize_plane_wave(a, lambda, b2, 1, tmp);
        // ONE source phase for both arrivals: that is what coherent means.
        const Real p = static_cast<Real>(2.0 * kPiD * rng.uniform01());
        const Complex rot(std::cos(p), std::sin(p));
        for (std::size_t i = 0; i < kElem; ++i) {
            g_snapshot[i] = (g_snapshot[i] + tmp[i] * static_cast<Real>(0.9)) * rot
                          + Complex(static_cast<Real>(0.05 * rng.normal()),
                                    static_cast<Real>(0.05 * rng.normal()));
        }
        covariance_accumulate(std::span<const Complex>(g_snapshot.data(), kElem), kElem, g_cov);
    }
    covariance_normalise(kElem, snaps, g_cov);

    constexpr std::size_t kAngles = 1801;
    const Real lo = deg2rad(static_cast<Real>(-20));
    const Real hi = deg2rad(static_cast<Real>(20));

    auto peaks_near = [&](std::span<const Real> p) {
        double top = 0;
        for (const Real v : p) top = std::max(top, static_cast<double>(v));
        std::size_t found = 0;
        for (std::size_t i = 1; i + 1 < p.size(); ++i) {
            if (p[i] > p[i - 1] && p[i] > p[i + 1] && static_cast<double>(p[i]) > top * 0.25) {
                ++found;
            }
        }
        return found;
    };

    // Unsmoothed MVDR on coherent sources.
    static std::array<Real, 4096> raw{};
    PT_CHECK(mvdr_power(a, lambda, g_cov, static_cast<Real>(0.001), lo, hi, g_mvdr_work,
                        std::span<Real>(raw.data(), kAngles)) == kAngles);
    const std::size_t raw_peaks = peaks_near(std::span<const Real>(raw.data(), kAngles));

    // Forward-backward smoothing over subarrays of 10 of the 16 elements.
    constexpr std::size_t kSub = 10;
    static std::array<Complex, kSub * kSub> smoothed{};
    static std::array<Complex, kSub * kSub + 2 * kSub> sub_work{};
    PT_CHECK(spatial_smooth(g_cov, kElem, kSub, smoothed));

    LineArray sub = a;
    sub.element_count = kSub;
    static std::array<Real, 4096> sm{};
    PT_CHECK(mvdr_power(sub, lambda, smoothed, static_cast<Real>(0.001), lo, hi, sub_work,
                        std::span<Real>(sm.data(), kAngles)) == kAngles);
    const std::size_t sm_peaks = peaks_near(std::span<const Real>(sm.data(), kAngles));

    std::printf("       two COHERENT sources at -6 and +6 deg\n");
    std::printf("       plain MVDR (16 elements)            : %zu peak(s)\n", raw_peaks);
    std::printf("       forward-backward smoothed (%zu-element subarrays): %zu peak(s)\n",
                kSub, sm_peaks);
    PT_CHECK(sm_peaks == 2);

    // The smoothed peaks land on the sources.
    double top = 0;
    for (std::size_t i = 0; i < kAngles; ++i) top = std::max(top, static_cast<double>(sm[i]));
    double found[2] = {0, 0};
    std::size_t k = 0;
    for (std::size_t i = 1; i + 1 < kAngles && k < 2; ++i) {
        if (sm[i] > sm[i - 1] && sm[i] > sm[i + 1] && static_cast<double>(sm[i]) > top * 0.25) {
            found[k++] = -20.0 + 40.0 * static_cast<double>(i) / static_cast<double>(kAngles - 1);
        }
    }
    std::printf("       smoothed peaks at %+.2f and %+.2f deg (truth -6, +6)\n",
                found[0], found[1]);
    PT_CHECK_NEAR(found[0], -6.0, 1.2);
    PT_CHECK_NEAR(found[1], 6.0, 1.2);

    // The cost, stated: smoothing spends aperture. A 10-element subarray has the
    // resolution of a 10-element array, not a 16-element one.
    std::printf("       cost: resolution falls from %.2f deg (16 elements) to %.2f (10)\n",
                static_cast<double>(rad2deg(conventional_resolution_limit_rad(a, lambda, 0))),
                static_cast<double>(rad2deg(conventional_resolution_limit_rad(sub, lambda, 0))));
    PT_CHECK(conventional_resolution_limit_rad(sub, lambda, 0)
           > conventional_resolution_limit_rad(a, lambda, 0));

    // Rejections.
    PT_CHECK(!spatial_smooth(g_cov, kElem, 0, smoothed));
    PT_CHECK(!spatial_smooth(g_cov, kElem, kElem + 1, smoothed));
    PT_CHECK(!spatial_smooth(g_cov, kElem, kSub, std::span<Complex>(smoothed.data(), 4)));
}

PT_TEST(spatial_smoothing_preserves_a_single_source) {
    // Smoothing must not move a bearing it had no business touching.
    const Real lambda = wavelength_m(10000, kC);
    const LineArray a = half_wave(kElem, lambda);
    const Real truth = deg2rad(static_cast<Real>(9));

    pt::Rng rng(303);
    covariance_clear(kElem, g_cov);
    for (std::size_t t = 0; t < 300; ++t) {
        synthesize_plane_wave(a, lambda, truth, 1, g_snapshot);
        const Real p = static_cast<Real>(2.0 * kPiD * rng.uniform01());
        const Complex rot(std::cos(p), std::sin(p));
        for (std::size_t i = 0; i < kElem; ++i) {
            g_snapshot[i] = g_snapshot[i] * rot
                          + Complex(static_cast<Real>(0.1 * rng.normal()),
                                    static_cast<Real>(0.1 * rng.normal()));
        }
        covariance_accumulate(std::span<const Complex>(g_snapshot.data(), kElem), kElem, g_cov);
    }
    covariance_normalise(kElem, 300, g_cov);

    constexpr std::size_t kSub = 12;
    static std::array<Complex, kSub * kSub> smoothed{};
    static std::array<Complex, kSub * kSub + 2 * kSub> sub_work{};
    PT_CHECK(spatial_smooth(g_cov, kElem, kSub, smoothed));

    LineArray sub = a;
    sub.element_count = kSub;
    constexpr std::size_t kAngles = 1801;
    const Real lo = deg2rad(static_cast<Real>(-30));
    const Real hi = deg2rad(static_cast<Real>(30));
    mvdr_power(sub, lambda, smoothed, static_cast<Real>(0.01), lo, hi, sub_work,
               std::span<Real>(g_power.data(), kAngles));
    const double est = static_cast<double>(rad2deg(estimate_bearing_rad(
        std::span<const Real>(g_power.data(), kAngles), lo, hi)));
    std::printf("       single source at 9.00 deg, smoothed MVDR reads %.3f\n", est);
    PT_CHECK_NEAR(est, 9.0, 0.4);
}
