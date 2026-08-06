// SPDX-License-Identifier: Apache-2.0
//
// A passive intercept receiver: two seconds of hydrophone data containing
// several active sonar pings, processed in a streaming overlap-save loop, with
// each detection reported as a Pulse Descriptor Word and scored against truth.
//
// The streaming pattern here is the part worth copying. The library never
// allocates and never holds state between calls, so the caller owns the sample
// history and slides it forward by exactly `bank.stride()` samples per block.
//
//   ./ping_intercept [output_dir]
#include "phantom/ping_analyzer.hpp"

#include <array>
#include <chrono>
#include <cstdio>
#include <cstdint>
#include <cmath>
#include <span>
#include <string>

using namespace phantom;

namespace {

constexpr Real        kFs       = 96000;
constexpr std::size_t kFft      = 8192;
constexpr Real        kDuration = 2;
constexpr std::size_t kTotal    = static_cast<std::size_t>(kDuration * kFs);

// Everything is statically allocated, exactly as it would be on the target.
std::array<Real, kTotal>  g_stream;
std::array<Real, kFft>    g_block;
PulseBank<8, kFft>        g_bank(kFs);
AnalyzerScratch<kFft>     g_scratch;

struct Truth {
    Real time_s;
    PulseType type;
    Real amplitude;
    const char* note;
};

PulseSpec make(PulseType t, Real f0, Real f1, Real dur) {
    PulseSpec s;
    s.type = t;
    s.f_start_hz = f0;
    s.f_end_hz = f1;
    s.duration_s = dur;
    return s;
}

// Deterministic noise, so the example prints the same thing on every machine.
struct Rng {
    std::uint64_t s = 0x9E3779B97F4A7C15ULL;
    double normal() noexcept {
        // Sum of twelve uniforms: close enough to Gaussian for a demo, and it
        // needs no transcendental functions.
        double acc = 0;
        for (int i = 0; i < 12; ++i) {
            s = s * 6364136223846793005ULL + 1442695040888963407ULL;
            acc += static_cast<double>((s >> 11) & 0xFFFFFF) / 16777216.0;
        }
        return acc - 6.0;
    }
};

void add_ping(const PulseSpec& spec, Real at_s, Real amplitude) {
    static std::array<Real, kFft> tmp{};
    const std::size_t n = render_real(spec, kFs, tmp);
    const auto at = static_cast<std::size_t>(at_s * kFs);
    for (std::size_t i = 0; i < n && at + i < kTotal; ++i) {
        g_stream[at + i] += amplitude * tmp[i];
    }
}

}  // namespace

int main(int argc, char** argv) {
    const std::string out_dir = (argc > 1) ? std::string(argv[1]) : std::string(".");

    std::printf("libphantom-sonar %s  --  passive ping intercept\n\n", PHANTOM_VERSION_STRING);

    // ---- 1. Waveform hypotheses -------------------------------------------
    // The bank must span what you expect to intercept. A waveform outside it
    // is not detected coherently; at best it shows up as a smeared partial
    // match on the nearest template.
    const PulseSpec cw       = make(PulseType::Cw,      12000, 12000, static_cast<Real>(0.020));
    const PulseSpec lfm_up   = make(PulseType::LfmUp,    8000, 20000, static_cast<Real>(0.020));
    const PulseSpec lfm_down = make(PulseType::LfmDown, 20000,  8000, static_cast<Real>(0.020));
    const PulseSpec hfm      = make(PulseType::Hfm,      8000, 20000, static_cast<Real>(0.020));

    g_bank.add(cw);
    g_bank.add(lfm_up);
    g_bank.add(lfm_down);
    g_bank.add(hfm);

    std::printf("Filter bank\n");
    for (std::size_t i = 0; i < g_bank.size(); ++i) {
        const PulseSpec& s = g_bank.view().templates[i].spec;
        std::printf("  [%zu] %-9s %6.0f - %6.0f Hz, %5.1f ms, TB = %6.0f\n",
                    i, pulse_type_name(s.type),
                    static_cast<double>(s.f_start_hz), static_cast<double>(s.f_end_hz),
                    static_cast<double>(s.duration_s) * 1e3,
                    static_cast<double>(s.time_bandwidth_product()));
    }
    std::printf("  transform %zu points, block advance %zu samples (%.1f ms)\n",
                kFft, g_bank.stride(),
                static_cast<double>(g_bank.stride()) / static_cast<double>(kFs) * 1e3);

    // ---- 2. Synthetic reception -------------------------------------------
    const Truth truths[] = {
        {static_cast<Real>(0.20), PulseType::LfmUp,   static_cast<Real>(1.00), "strong direct path"},
        {static_cast<Real>(0.55), PulseType::Cw,      static_cast<Real>(0.60), "CW searchlight"},
        {static_cast<Real>(0.95), PulseType::Hfm,     static_cast<Real>(0.35), "HFM, Doppler tolerant"},
        {static_cast<Real>(1.30), PulseType::LfmDown, static_cast<Real>(0.25), "downsweep, weak"},
        {static_cast<Real>(1.70), PulseType::LfmUp,   static_cast<Real>(0.15), "distant repeat"},
    };
    const PulseSpec* specs[] = {&lfm_up, &cw, &hfm, &lfm_down, &lfm_up};

    Rng rng;
    const double noise_sigma = 0.5;
    for (Real& v : g_stream) v = static_cast<Real>(noise_sigma * rng.normal());
    for (std::size_t i = 0; i < 5; ++i) add_ping(*specs[i], truths[i].time_s, truths[i].amplitude);

    std::printf("\nReceived stream\n");
    std::printf("  %.1f s at %.0f kHz, noise sigma %.2f, %zu pings\n",
                static_cast<double>(kDuration), static_cast<double>(kFs) / 1000.0,
                noise_sigma, sizeof(truths) / sizeof(truths[0]));

    // ---- 3. Detector setup -------------------------------------------------
    DetectorConfig cfg;
    // The guard band must clear the widest response in the bank. The CW does
    // not pulse-compress, so its correlation is as wide as the pulse itself;
    // a smaller guard makes the CW mask itself and vanish silently.
    cfg.cfar_guard = suggested_cfar_guard(g_bank.view());
    cfg.cfar_train = 256;
    cfg.threshold_alpha = cfar_alpha(2 * cfg.cfar_train, static_cast<Real>(1e-6));
    cfg.dead_time_s = suggested_dead_time_s(g_bank.view());

    std::printf("\nDetector\n");
    std::printf("  CA-CFAR, guard %zu cells, %zu training cells either side\n",
                cfg.cfar_guard, cfg.cfar_train);
    std::printf("  design Pfa 1e-6 -> threshold %.1f dB over the local noise\n",
                10.0 * std::log10(static_cast<double>(cfg.threshold_alpha)));
    std::printf("  dead time %.1f ms (the longest replica; a shorter one splits the\n"
                "  uncompressed CW response into repeated detections)\n",
                static_cast<double>(cfg.dead_time_s) * 1e3);

    // ---- 4. Streaming loop -------------------------------------------------
    std::array<PulseDescriptor, 16> pdw{};
    const std::size_t stride = g_bank.stride();

    std::printf("\nDetections\n");
    std::printf("  %-10s %-9s %8s %9s %10s %9s\n",
                "ToA (s)", "type", "amp", "SNR dB", "fc (Hz)", "BW (Hz)");

    std::size_t reported = 0;
    std::size_t matched = 0;
    double worst_toa_err_us = 0.0;
    bool truth_found[5] = {false, false, false, false, false};

    const auto t0 = std::chrono::steady_clock::now();
    std::size_t blocks = 0;

    for (std::size_t start = 0; start + kFft <= kTotal; start += stride) {
        // Overlap-save: hand the analyser a full transform's worth of samples,
        // but only lags [0, stride) are alias-free, so the window advances by
        // `stride` and the trailing (kFft - stride) samples are re-presented.
        for (std::size_t i = 0; i < kFft; ++i) g_block[i] = g_stream[start + i];
        const Real block_time = static_cast<Real>(start) / kFs;

        const std::size_t n = analyze_block(g_bank.view(), cfg, g_block, block_time,
                                            g_scratch.view(), pdw);
        ++blocks;

        for (std::size_t i = 0; i < n; ++i) {
            ++reported;
            std::printf("  %-10.5f %-9s %8.3f %9.1f %10.0f %9.0f",
                        static_cast<double>(pdw[i].toa_s),
                        pulse_type_name(pdw[i].type),
                        static_cast<double>(pdw[i].amplitude),
                        static_cast<double>(pdw[i].snr_db),
                        static_cast<double>(pdw[i].centre_freq_hz),
                        static_cast<double>(pdw[i].bandwidth_hz));

            // Score against truth. The gate is waveform-dependent on purpose:
            // range resolution is 1/B, so a 12 kHz chirp localises an arrival
            // to tens of microseconds while a CW -- zero bandwidth -- can only
            // place it to within a fraction of its own length. That is not a
            // detector defect, it is why chirps are transmitted.
            bool hit = false;
            for (std::size_t k = 0; k < 5; ++k) {
                const double err_s = static_cast<double>(pdw[i].toa_s)
                                   - static_cast<double>(truths[k].time_s);
                const double gate = (truths[k].type == PulseType::Cw)
                                  ? 0.5 * static_cast<double>(specs[k]->duration_s)
                                  : 1e-3;
                if (std::fabs(err_s) < gate && pdw[i].type == truths[k].type) {
                    if (!truth_found[k]) { truth_found[k] = true; ++matched; }
                    if (std::fabs(err_s) * 1e6 > worst_toa_err_us) {
                        worst_toa_err_us = std::fabs(err_s) * 1e6;
                    }
                    std::printf("   <- %s, ToA error %+.1f us\n", truths[k].note, err_s * 1e6);
                    hit = true;
                    break;
                }
            }
            if (!hit) std::printf("   <- unmatched\n");
        }
    }
    const auto t1 = std::chrono::steady_clock::now();

    const auto elapsed_us = static_cast<double>(
        std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count());

    // ---- 5. Summary --------------------------------------------------------
    std::printf("\nSummary\n");
    std::printf("  pings present        : 5\n");
    std::printf("  correctly detected   : %zu\n", matched);
    std::printf("  detections reported  : %zu\n", reported);
    std::printf("  worst ToA error      : %.1f us  (%.1f cm of one-way range)\n",
                worst_toa_err_us, worst_toa_err_us * 1e-6 * 1500.0 * 100.0);
    std::printf("                         chirps localise to ~1/B = %.0f us; the CW,\n"
                "                         having no bandwidth, cannot do better than a\n"
                "                         fraction of its own %.0f ms length\n",
                1e6 / 12000.0, static_cast<double>(cw.duration_s) * 1e3);
    for (std::size_t k = 0; k < 5; ++k) {
        if (!truth_found[k]) {
            std::printf("  MISSED               : %s at %.2f s (amplitude %.2f)\n",
                        pulse_type_name(truths[k].type),
                        static_cast<double>(truths[k].time_s),
                        static_cast<double>(truths[k].amplitude));
        }
    }
    if (reported > matched) {
        std::printf("  extra detections     : %zu -- cross-template ghosts and false\n"
                    "                         alarms; the LFM and HFM share a band, so an\n"
                    "                         arrival of one partially matches the other.\n",
                    reported - matched);
    }

    std::printf("\nTiming\n");
    std::printf("  %zu blocks, %.1f ms of audio\n", blocks, static_cast<double>(kDuration) * 1e3);
    std::printf("  processing           : %.1f ms wall, %.2f ms per block\n",
                elapsed_us / 1000.0, elapsed_us / 1000.0 / static_cast<double>(blocks));
    std::printf("  real-time factor     : %.0f x on one core\n",
                static_cast<double>(kDuration) * 1e6 / elapsed_us);
    std::printf("  detection latency    : %.1f ms worst case (one block advance)\n",
                static_cast<double>(stride) / static_cast<double>(kFs) * 1e3);

    // ---- 6. CSV for plotting ----------------------------------------------
    const std::string path = out_dir + "/intercept.csv";
    if (std::FILE* f = std::fopen(path.c_str(), "w")) {
        std::fprintf(f, "time_s,sample\n");
        for (std::size_t i = 0; i < kTotal; i += 4) {
            std::fprintf(f, "%.6f,%.5f\n", static_cast<double>(i) / static_cast<double>(kFs),
                         static_cast<double>(g_stream[i]));
        }
        std::fclose(f);
        std::printf("\nWrote %s\n", path.c_str());
    }
    return 0;
}
