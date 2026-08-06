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
#include "phantom/ping_analyzer.hpp"
#include "phantom/profile.hpp"
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
        {static_cast<Real>(8),  static_cast<Real>(-2),  static_cast<Real>(3), 1},
        {static_cast<Real>(20), static_cast<Real>(-5),  static_cast<Real>(0), static_cast<Real>(1.4)},
        {static_cast<Real>(32), static_cast<Real>(-8),  static_cast<Real>(-4), 1},
        {static_cast<Real>(44), static_cast<Real>(-11), static_cast<Real>(6), 1},
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

    // ---- 6. The claim this library will not make ---------------------------
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
