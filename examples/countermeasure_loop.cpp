// SPDX-License-Identifier: Apache-2.0
//
// The whole library in one loop: an active sonar pings, this vehicle
// intercepts and classifies the ping, and synthesises a ghost swarm in reply.
// The transmitted swarm is then fed back through the analyser -- standing in
// for the sonar operator -- so the result can be scored rather than admired.
//
// The ocean model supplies the sound speed, so the delays are computed with the
// speed that actually applies at depth rather than a textbook 1500 m/s.
//
//   ./countermeasure_loop
#include "phantom/echo_synth.hpp"
#include "phantom/boundary.hpp"
#include "phantom/eigenray.hpp"
#include "phantom/ping_analyzer.hpp"
#include "phantom/profile.hpp"
#include "phantom/reverberation.hpp"
#include "phantom/sound_speed.hpp"

#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>

using namespace phantom;

namespace {

constexpr Real        kFs  = 96000;
constexpr std::size_t kFft = 8192;

// Static storage throughout. A Doppler bank is megabytes; putting one on the
// stack is the first mistake to make with this API.
PulseBank<64, kFft>       g_bank(kFs);
AnalyzerScratch<kFft>     g_scratch;
std::array<Real, kFft>    g_rx;
std::array<Real, kFft>    g_tx;
SoundSpeedProfile<64>     g_profile;

struct Rng {
    std::uint64_t s = 0xD1B54A32D192ED03ULL;
    double normal() noexcept {
        double acc = 0;
        for (int i = 0; i < 12; ++i) {
            s = s * 6364136223846793005ULL + 1442695040888963407ULL;
            acc += static_cast<double>((s >> 11) & 0xFFFFFF) / 16777216.0;
        }
        return acc - 6.0;
    }
};

PulseSpec make(PulseType t, Real f0, Real f1, Real dur) {
    PulseSpec s;
    s.type = t;
    s.f_start_hz = f0;
    s.f_end_hz = f1;
    s.duration_s = dur;
    return s;
}

void place(std::span<const Real> sig, std::size_t at, std::span<Real> dst, Real gain = 1) {
    for (std::size_t i = 0; i < sig.size() && at + i < dst.size(); ++i) {
        dst[at + i] += gain * sig[i];
    }
}

}  // namespace

int main() {
    std::printf("libphantom-sonar %s  --  intercept and reply\n\n", PHANTOM_VERSION_STRING);

    // ---- 1. Ocean: the sound speed that sets every delay -------------------
    // A summer thermocline, and the vehicle sitting at 120 m.
    const Real vehicle_depth_m = 120;
    fill_profile(g_profile, 0, 300, 31, [](Real z) {
        return sound_speed::mackenzie(static_cast<Real>(18) - static_cast<Real>(0.03) * z, 35, z);
    });
    const Real c = speed_at(g_profile.view(), vehicle_depth_m);

    std::printf("Ocean\n");
    std::printf("  vehicle depth      : %.0f m\n", static_cast<double>(vehicle_depth_m));
    std::printf("  local sound speed  : %.2f m/s  (Mackenzie, %.1f C at depth)\n",
                static_cast<double>(c),
                static_cast<double>(18 - 0.03 * static_cast<double>(vehicle_depth_m)));
    std::printf("  a 1500 m/s assumption would misplace a 100 m ghost by %.2f m\n",
                100.0 * (static_cast<double>(c) / 1500.0 - 1.0));

    // ---- 2. The bank we listen with ----------------------------------------
    const PulseSpec cw  = make(PulseType::Cw,   12000, 12000, static_cast<Real>(0.02));
    const PulseSpec lfm = make(PulseType::LfmUp, 8000, 20000, static_cast<Real>(0.02));
    const PulseSpec hfm = make(PulseType::Hfm,   8000, 20000, static_cast<Real>(0.02));

    g_bank.add(cw);
    g_bank.add(hfm);
    const std::size_t lfm_bins = g_bank.add_doppler_bank(lfm, -20, 20, static_cast<Real>(0.5));

    std::printf("\nFilter bank\n");
    std::printf("  %zu templates: CW, HFM, and %zu LFM Doppler bins over +/- 20 m/s\n",
                g_bank.size(), lfm_bins);
    std::printf("  bin spacing        : %.2f m/s\n",
                static_cast<double>(g_bank.doppler_bin_spacing_mps(lfm, -20, 20,
                                                                   static_cast<Real>(0.5))));
    std::printf("  an HFM would need  : %zu bin(s) for the same span -- the design\n"
                "                       argument for hyperbolic sweeps in one number\n",
                g_bank.doppler_bins_required(hfm, -20, 20, static_cast<Real>(0.5)));
    std::printf("  memory             : %.1f MB of replica spectra\n",
                64.0 * static_cast<double>(kFft) * static_cast<double>(sizeof(Complex)) / 1048576.0);

    DetectorConfig cfg;
    cfg.cfar_guard = suggested_cfar_guard(g_bank.view());
    cfg.cfar_train = 256;
    cfg.threshold_alpha = cfar_alpha(512, static_cast<Real>(1e-6));
    cfg.dead_time_s = suggested_dead_time_s(g_bank.view());

    // ---- 3. Intercept -------------------------------------------------------
    // An LFM ping from a sonar closing at 8 m/s, arriving 5 ms into the block.
    const Real closing_mps = 8;
    Rng rng;
    for (Real& v : g_rx) v = static_cast<Real>(0.4 * rng.normal());

    static std::array<Real, kFft> tmp{};
    const std::size_t n = render_real_doppler(lfm, kFs, closing_mps / c, tmp);
    const auto ping_at = static_cast<std::size_t>(static_cast<Real>(0.005) * kFs);
    place(std::span<const Real>(tmp.data(), n), ping_at, g_rx, static_cast<Real>(0.9));

    std::array<PulseDescriptor, 16> pdw{};
    const auto t0 = std::chrono::steady_clock::now();
    const std::size_t got = analyze_block(g_bank.view(), cfg, g_rx, 0, g_scratch.view(), pdw);
    const auto t1 = std::chrono::steady_clock::now();

    std::printf("\nIntercept\n");
    std::printf("  transmitted        : LFM 8-20 kHz, 20 ms, closing at %.1f m/s\n",
                static_cast<double>(closing_mps));
    if (got == 0) {
        std::printf("  nothing detected -- nothing to reply to\n");
        return 1;
    }
    const PulseDescriptor& ping = pdw[0];
    std::printf("  detected           : %s, %.0f-%.0f Hz, %.1f ms\n",
                pulse_type_name(ping.type),
                static_cast<double>(ping.centre_freq_hz - ping.bandwidth_hz / 2),
                static_cast<double>(ping.centre_freq_hz + ping.bandwidth_hz / 2),
                static_cast<double>(ping.duration_s) * 1e3);
    const double toa_truth_ms = static_cast<double>(ping_at) / static_cast<double>(kFs) * 1e3;
    const double toa_err_ms = static_cast<double>(ping.toa_s) * 1e3 - toa_truth_ms;
    const double bin_spacing = static_cast<double>(
        g_bank.doppler_bin_spacing_mps(lfm, -20, 20, static_cast<Real>(0.5)));
    std::printf("  radial velocity    : %.1f m/s%s (truth %.1f, bins %.1f m/s apart\n"
                "                       so a half-bin error is expected)\n",
                static_cast<double>(ping.radial_velocity_mps),
                ping.doppler_resolved ? "" : " [unresolved]",
                static_cast<double>(closing_mps), bin_spacing);
    std::printf("  arrival            : %.4f ms (truth %.4f, error %+.4f)\n",
                static_cast<double>(ping.toa_s) * 1e3, toa_truth_ms, toa_err_ms);

    // That arrival error is not noise. The winning Doppler bin is off by
    // v_err, and an LFM converts residual Doppler into a range bias by
    // dt = -(v_err/c) * f_end / mu -- the wideband coupling derived in
    // docs/math_spec.md 6.2. Predicting it here is the cheapest possible check
    // that the two halves of the library agree with each other.
    {
        const double v_err = static_cast<double>(ping.radial_velocity_mps)
                           - static_cast<double>(closing_mps);
        const double mu = static_cast<double>(lfm.chirp_rate_hz_s());
        const double f_end = static_cast<double>(lfm.f_end_hz);
        const double predicted_ms = (v_err / static_cast<double>(c)) * f_end / mu * 1e3;
        std::printf("                       predicted from the %+.1f m/s bin error:\n"
                    "                       dt = -(v/c) f_end/mu = %+.4f ms\n",
                    v_err, predicted_ms);
    }
    std::printf("  SNR                : %.1f dB\n", static_cast<double>(ping.snr_db));
    std::printf("  analysis time      : %.2f ms for a %zu-template bank\n",
                static_cast<double>(
                    std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count()) / 1e3,
                g_bank.size());

    // ---- 4. Reply -----------------------------------------------------------
    // Four ghosts: one near the true position and three spread beyond it, with
    // decreasing target strength so the set looks like a dispersed formation.
    const EchoSpec swarm[] = {
        {.range_offset_m = 8,  .target_strength_db = -2,  .radial_velocity_mps = 3},
        {.range_offset_m = 20, .target_strength_db = -5,  .length_scale = static_cast<Real>(1.4)},
        {.range_offset_m = 32, .target_strength_db = -8,  .radial_velocity_mps = -4},
        {.range_offset_m = 44, .target_strength_db = -11, .radial_velocity_mps = 6},
    };
    const std::span<const EchoSpec> echoes(swarm, 4);

    const std::size_t need = swarm_length(ping, echoes, kFs, c);
    const std::size_t written = synthesize_swarm(ping, echoes, kFs, c, g_tx);

    std::printf("\nReply\n");
    std::printf("  ghosts             : %zu, spanning %.1f ms of transmission\n",
                echoes.size(), static_cast<double>(need) / static_cast<double>(kFs) * 1e3);
    std::printf("  %10s %10s %10s %10s\n", "range (m)", "TS (dB)", "v (m/s)", "delay (ms)");
    for (const EchoSpec& e : echoes) {
        std::printf("  %10.0f %10.0f %10.0f %10.2f\n",
                    static_cast<double>(e.range_offset_m),
                    static_cast<double>(e.target_strength_db),
                    static_cast<double>(e.radial_velocity_mps),
                    static_cast<double>(echo_delay_s(e.range_offset_m, c)) * 1e3);
    }
    if (written == 0) {
        std::printf("  transmit buffer too short\n");
        return 1;
    }

    // ---- 5. What the sonar operator sees ------------------------------------
    for (Real& v : g_rx) v = static_cast<Real>(0.4 * rng.normal());
    place(std::span<const Real>(g_tx.data(), written), 200, g_rx);

    DetectorConfig rx_cfg = cfg;
    rx_cfg.dead_time_s = static_cast<Real>(0.004);
    const std::size_t seen = analyze_block(g_bank.view(), rx_cfg, g_rx, 0, g_scratch.view(), pdw);

    std::printf("\nAs received back\n");
    std::printf("  %12s %-9s %9s %9s %10s\n",
                "ToA (ms)", "type", "amp", "SNR dB", "v (m/s)");
    std::size_t matched = 0;
    for (std::size_t i = 0; i < seen; ++i) {
        const double toa_ms = static_cast<double>(pdw[i].toa_s) * 1e3;
        const double rel_ms = toa_ms - 200.0 / static_cast<double>(kFs) * 1e3;
        std::printf("  %12.3f %-9s %9.3f %9.1f %10.1f",
                    toa_ms, pulse_type_name(pdw[i].type),
                    static_cast<double>(pdw[i].amplitude),
                    static_cast<double>(pdw[i].snr_db),
                    static_cast<double>(pdw[i].radial_velocity_mps));
        bool hit = false;
        for (const EchoSpec& e : echoes) {
            const double want = static_cast<double>(echo_delay_s(e.range_offset_m, c)) * 1e3;
            if (std::fabs(rel_ms - want) < 0.2) {
                std::printf("   <- ghost at %.0f m\n", static_cast<double>(e.range_offset_m));
                ++matched;
                hit = true;
                break;
            }
        }
        if (!hit) std::printf("   <- unresolved\n");
    }

    std::printf("\nSummary\n");
    std::printf("  ghosts transmitted : %zu\n", echoes.size());
    std::printf("  ghosts resolved    : %zu\n", matched);
    std::printf("  total detections   : %zu\n", seen);
    if (matched < echoes.size()) {
        std::printf("\n  The 20 m ghost carries length_scale = 1.4, so it is a 28 ms pulse\n"
                    "  and the bank holds no 28 ms template. It is detected, but by a\n"
                    "  mismatched replica, so its reported time and type are both wrong.\n"
                    "  That cuts both ways: an extended-target echo is more convincing to\n"
                    "  an operator AND harder for any bank -- including this one -- to\n"
                    "  pin down. A bank only resolves what it spans.\n");
    }

    // ---- 6. The multipath the ocean itself imposes -------------------------
    // Everything above placed ghosts at chosen delays. The ray tracer says what
    // a single reflector at that range ACTUALLY returns: several arrivals, at
    // times and levels neither of us picked.
    static std::array<RayPoint, 8192> ray_scratch{};
    static std::array<Eigenray, 24> paths{};
    TraceConfig trace_cfg;
    trace_cfg.max_range_m = 20000;
    trace_cfg.max_time_s = 60;
    trace_cfg.bottom_depth_m = 300;

    EigenraySearch es;
    es.angle_min_rad = deg2rad(static_cast<Real>(-10));
    es.angle_max_rad = deg2rad(static_cast<Real>(10));
    es.fan_count = 1201;

    const Real sonar_range_m = 2500;
    const std::size_t np = find_eigenrays(g_profile.view(), vehicle_depth_m,
                                          static_cast<Real>(40), sonar_range_m,
                                          trace_cfg, es, ray_scratch, paths);

    std::printf("\nMultipath from the ray tracer\n");
    std::printf("  sonar at 40 m, this vehicle at %.0f m, %.1f km apart\n",
                static_cast<double>(vehicle_depth_m),
                static_cast<double>(sonar_range_m) / 1000.0);
    std::printf("  %zu eigenrays found\n", np);
    // Boundary losses. Before v0.5 every reflection was perfect, which made a
    // four-bounce path look as loud as the direct one -- the largest remaining
    // overstatement in the transmission loss.
    BoundaryModel bounds;
    bounds.surface.rms_wave_height_m = wind_to_rms_wave_height_m(8);   // 8 m/s wind
    const Real c_surface = speed_at(g_profile.view(), 0);
    const Real c_bottom = speed_at(g_profile.view(), 300);
    const Real crit = bottom_critical_angle_rad(bounds.bottom, c_bottom);

    std::printf("  sediment %.0f m/s, density x%.1f, %.1f dB/wavelength\n",
                static_cast<double>(bounds.bottom.sound_speed_mps),
                static_cast<double>(bounds.bottom.density_ratio),
                static_cast<double>(bounds.bottom.attenuation_db_per_wavelength));
    std::printf("  critical grazing angle %.2f deg -- below it the bottom traps,\n"
                "  above it energy leaks into the sediment\n",
                static_cast<double>(rad2deg(crit)));
    std::printf("  sea state: 8 m/s wind, %.2f m RMS wave height\n",
                static_cast<double>(bounds.surface.rms_wave_height_m));

    std::printf("  %10s %10s %8s %8s %10s %12s %10s\n",
                "launch", "t (ms)", "srf", "btm", "graze", "TL no bnd", "TL total");
    for (std::size_t i = 0; i < np; ++i) {
        const double plain = static_cast<double>(
            transmission_loss_db(paths[i], sonar_range_m, ping.centre_freq_hz, c));
        const double total = static_cast<double>(total_transmission_loss_db(
            paths[i], sonar_range_m, ping.centre_freq_hz, c, bounds, c_surface, c_bottom));
        const double graze = static_cast<double>(rad2deg(static_cast<Real>(std::acos(
            std::min(1.0, static_cast<double>(paths[i].snell_invariant)
                        * static_cast<double>(c_bottom))))));
        std::printf("  %9.3f%s %10.3f %8u %8u %9.2f%s %12.2f %10.2f%s\n",
                    static_cast<double>(rad2deg(paths[i].launch_angle_rad)), "d",
                    static_cast<double>(paths[i].travel_time_s) * 1e3,
                    paths[i].surface_bounces, paths[i].bottom_bounces,
                    graze, "d", plain, total,
                    paths[i].near_caustic ? "  [caustic]" : "");
    }

    if (np > 0) {
        static std::array<EchoSpec, 24> mp{};
        const std::size_t nm = echoes_from_eigenrays(
            std::span<const Eigenray>(paths.data(), np), sonar_range_m,
            ping.centre_freq_hz, c, static_cast<Real>(-6), mp);
        double spread_ms = 0;
        for (std::size_t i = 0; i < nm; ++i) {
            spread_ms = std::max(spread_ms, static_cast<double>(mp[i].extra_delay_s) * 1e3);
        }
        std::printf("  -> %zu echoes spanning %.2f ms two-way\n", nm, spread_ms);
        std::printf("     The transmitted pulse is %.1f ms long, so these paths %s.\n",
                    static_cast<double>(ping.duration_s) * 1e3,
                    (spread_ms > static_cast<double>(ping.duration_s) * 1e3)
                        ? "resolve as separate arrivals"
                        : "smear one arrival rather than resolving");
        std::printf("     Absorption at %.0f kHz costs %.2f dB/km, so the longest path\n"
                    "     pays %.1f dB more than the shortest on distance alone.\n",
                    static_cast<double>(ping.centre_freq_hz) / 1000.0,
                    static_cast<double>(thorp_absorption_db_per_km(ping.centre_freq_hz)),
                    static_cast<double>(thorp_absorption_db_per_km(ping.centre_freq_hz))
                        * static_cast<double>(paths[np - 1].path_length_m
                                              - paths[0].path_length_m) / 1000.0);
    }

    // ---- 7. What the sonar is actually fighting ----------------------------
    // The analyser assumes white noise. A real active sonar at these ranges is
    // reverberation-limited, and that changes which knobs help.
    {
        const Real beam = static_cast<Real>(0.15);          // ~8.6 degrees
        const Real ts = 10;                                 // target strength
        const Real graze = deg2rad(static_cast<Real>(12));
        const Real ss_bottom = lambert_bottom_scattering_db(static_cast<Real>(-27), graze);
        const Real ss_surface = chapman_harris_surface_scattering_db(
            mps_to_knots(8), ping.centre_freq_hz, graze);

        std::printf("\nReverberation\n");
        std::printf("  bottom scattering (Lambert, 12 deg)  : %.1f dB\n",
                    static_cast<double>(ss_bottom));
        std::printf("  surface scattering (8 m/s, 12 deg)   : %.1f dB\n",
                    static_cast<double>(ss_surface));

        // The uncompressed pulse against its compressed equivalent. A chirp's
        // range resolution is 1/B, not its length, so the ensonified area --
        // and with it the reverberation -- shrinks by the time-bandwidth
        // product.
        const Real tau_raw = ping.duration_s;
        const Real tau_eff = (ping.bandwidth_hz > 0)
                           ? static_cast<Real>(1) / ping.bandwidth_hz : tau_raw;
        const Real area_raw = ensonified_area_m2(sonar_range_m, beam, tau_raw, c);
        const Real area_eff = ensonified_area_m2(sonar_range_m, beam, tau_eff, c);

        const double er_raw = static_cast<double>(
            echo_to_reverberation_ratio_db(ts, ss_bottom, area_raw));
        const double er_eff = static_cast<double>(
            echo_to_reverberation_ratio_db(ts, ss_bottom, area_eff));

        std::printf("  ensonified area at %.1f km          : %.0f m^2 raw, %.0f m^2 compressed\n",
                    static_cast<double>(sonar_range_m) / 1000.0,
                    static_cast<double>(area_raw), static_cast<double>(area_eff));
        std::printf("  echo-to-reverberation ratio          : %+.1f dB raw, %+.1f dB compressed\n",
                    er_raw, er_eff);
        std::printf("  pulse compression is worth            %+.1f dB here\n", er_eff - er_raw);

        std::printf("\n  Note what does NOT appear in that ratio: source level and\n");
        std::printf("  transmission loss. EL - RL = TS - S_s - 10log10(A), so they cancel.\n");
        std::printf("  Against reverberation a louder transmitter buys nothing at all --\n");
        std::printf("  only a shorter effective pulse or a narrower beam does. That is the\n");
        std::printf("  design argument for the chirps this analyser is built around.\n");

        const Real crossover = reverberation_limited_range_m(
            200, ss_bottom, beam, tau_eff, c, 60);
        std::printf("\n  Reverberation meets a 60 dB ambient at %.0f m: inside that the\n"
                    "  geometry is reverberation-limited, outside it noise-limited.\n",
                    static_cast<double>(crossover));
    }

    // ---- 8. The claim this library will not make ---------------------------
    const Real fc = ping.centre_freq_hz;
    const Real tau20 = max_timing_error_s(static_cast<Real>(-20), fc);
    const Real tau10 = max_timing_error_s(static_cast<Real>(-10), fc);
    std::printf("\nOn cancelling the echo instead of faking it\n");
    std::printf("  At %.0f kHz, anti-phase cancellation needs the timing held to\n",
                static_cast<double>(fc) / 1000.0);
    std::printf("    %.2f us (%.2f mm of path) for 10 dB\n",
                static_cast<double>(tau10) * 1e6, static_cast<double>(tau10 * c) * 1e3);
    std::printf("    %.2f us (%.2f mm of path) for 20 dB\n",
                static_cast<double>(tau20) * 1e6, static_cast<double>(tau20 * c) * 1e3);
    std::printf("  Past a quarter period the canceller adds %.1f dB instead of removing\n"
                "  any, and hull scattering is distributed over many wavelengths, so one\n"
                "  projector cannot match phase at more than one bearing. Generating\n"
                "  false targets is the achievable countermeasure; cancelling the real\n"
                "  one is not. See docs/validation.md.\n",
                static_cast<double>(null_gain_db(ping.bandwidth_hz, fc,
                                                 static_cast<Real>(1) / (4 * fc))));
    return 0;
}
