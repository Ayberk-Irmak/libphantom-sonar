// SPDX-License-Identifier: Apache-2.0
// libphantom-sonar — active sonar ping detection and characterisation.
//
// This is the stage the original project specification skipped. An echo
// synthesiser needs to know what arrived, when, and what shape it was; without
// that it has nothing to respond to. Detecting and characterising the incoming
// ping is also the harder half of the problem, and the one with the binding
// latency budget.
//
// Structure: a bank of matched filters, one per waveform hypothesis, run over
// the same block. Taking the largest normalised peak across the bank is a
// generalised likelihood ratio test over a discrete hypothesis set -- the
// standard estimator, and the reason the bank must span the waveforms you
// expect rather than just one.
//
// Detection is cell-averaging CFAR on the compressed output, so the threshold
// tracks the local noise level instead of being a fixed number that is wrong
// the moment the sea state changes.
//
// Streaming: the caller owns the sample history. Each block is `fft_size`
// samples; after processing, advance by `stride()` and keep the last
// `max_replica_length() - 1` samples. That is overlap-save, and it is why lags
// beyond the stride are not reported: they would be circular-correlation
// aliases.
#ifndef PHANTOM_PING_ANALYZER_HPP
#define PHANTOM_PING_ANALYZER_HPP

#include "phantom/fft.hpp"
#include "phantom/matched_filter.hpp"
#include "phantom/types.hpp"
#include "phantom/waveform.hpp"

#include <array>
#include <span>

namespace phantom {

// What the analyser reports about one detected pulse: a Pulse Descriptor Word.
struct PulseDescriptor {
    Real toa_s = 0;             // arrival time of the pulse leading edge
    Real snr_db = 0;            // post-compression SNR over the local noise
    Real amplitude = 0;         // estimated received amplitude
    Real peak_magnitude = 0;    // raw matched filter peak

    PulseType type = PulseType::Unknown;
    Real centre_freq_hz = 0;
    Real bandwidth_hz = 0;
    Real duration_s = 0;
    Real chirp_rate_hz_s = 0;

    // Radial velocity of the transmitter relative to this receiver, from the
    // Doppler bin that won. Zero when the bank has no Doppler coverage --
    // which is not the same as "the target is stationary", so check
    // `doppler_resolved` before believing it.
    Real radial_velocity_mps = 0;
    bool doppler_resolved = false;

    std::size_t template_index = 0;
};

struct PulseTemplate {
    PulseSpec spec;
    MatchedFilter filter;
    // Time-scale offset this replica was rendered at: delta = v/c, positive
    // for a closing geometry.
    Real doppler = 0;
    bool doppler_bin = false;   // part of a deliberate Doppler sweep
};

struct AnalyzerView {
    std::span<const PulseTemplate> templates;
    FftView fft;
    Real sample_rate_hz = 0;
    std::size_t max_replica_length = 0;
    Real sound_speed_mps = 1500;

    [[nodiscard]] bool valid() const noexcept;
    // Usable lags per block, and the block advance. Set by the LONGEST replica
    // in the bank: a lag is only reported where every template is alias-free.
    [[nodiscard]] std::size_t stride() const noexcept;
};

struct DetectorConfig {
    // Linear multiplier on the CFAR noise-power estimate. Use cfar_alpha() to
    // derive it from a target false-alarm probability.
    Real threshold_alpha = 20;
    // Cells either side excluded from the noise estimate. THIS MUST EXCEED THE
    // WIDTH OF THE RESPONSE YOU ARE TRYING TO DETECT -- see
    // suggested_cfar_guard(), and read its note before choosing a value.
    std::size_t cfar_guard = 16;
    std::size_t cfar_train = 128;    // training cells either side
    Real dead_time_s = 0;            // minimum spacing between reported pulses
};

// CA-CFAR scaling for a target false-alarm probability, assuming exponentially
// distributed output power (Rayleigh magnitude):  alpha = N (Pfa^(-1/N) - 1).
[[nodiscard]] Real cfar_alpha(std::size_t training_cells, Real pfa) noexcept;

// Smallest guard band that keeps a bank's own responses out of its noise
// estimate: the length of the longest replica.
//
// This is not a refinement, it is a correctness requirement, and getting it
// wrong fails silently in the worst way. Cell-averaging CFAR estimates the
// noise from cells either side of the one under test. If the target's response
// is wider than the guard band, the target leaks into its own training cells,
// the threshold rises with the signal, and the ratio never crosses -- the
// stronger the pulse, the higher the threshold it has to beat. The detector
// reports nothing and looks like it is working.
//
// A chirp compresses to about fs/B cells, so a small guard suffices. A CW does
// not compress at all: its correlation is a triangle 2L wide. A bank holding
// both must be guarded for the CW, or it will detect every waveform except the
// simplest one.
[[nodiscard]] std::size_t suggested_cfar_guard(const AnalyzerView& bank) noexcept;

// Smallest dead time that stops one arrival being reported several times: the
// longest replica's duration. Same reasoning as suggested_cfar_guard(). A
// compressed chirp needs only a few cells, but a CW's response is as wide as
// the pulse, so a short dead time chops it into a burst of detections spaced
// exactly one dead time apart -- a distinctive and easily misread artefact.
[[nodiscard]] Real suggested_dead_time_s(const AnalyzerView& bank) noexcept;

// Caller-owned working buffers. Every span must be at least `fft_size` long
// (`prefix` one longer). Nothing is allocated inside the analyser.
struct AnalyzerWorkspace {
    std::span<Complex> fft_scratch;
    std::span<Complex> correlation;
    std::span<Real> best_power;             // max over the bank, per lag
    std::span<std::uint16_t> best_template; // which template won that lag
    std::span<Real> prefix;                 // running sums, for O(1) CFAR

    [[nodiscard]] bool valid(std::size_t fft_size) const noexcept;
};

// Analyses one block of `fft_size` real samples. `block_start_time_s` is the
// timestamp of `block[0]`. Returns the number of descriptors written.
std::size_t analyze_block(const AnalyzerView& bank,
                          const DetectorConfig& cfg,
                          std::span<const Real> block,
                          Real block_start_time_s,
                          const AnalyzerWorkspace& work,
                          std::span<PulseDescriptor> out) noexcept;

// Owning bank.
//
// SIZE THIS OBJECT DELIBERATELY, AND DO NOT PUT IT ON THE STACK. The spectra
// alone cost MaxTemplates * FftSize * sizeof(Complex): 2 MB at 16 templates and
// an 8192-point transform, 8.4 MB at 64. A Doppler bank reaches those counts
// easily -- an LFM over +/-20 m/s needs dozens of bins -- so a bank declared as
// a local variable will overflow a default thread stack before it ever runs.
// Declare it static, or as a member of something that lives in static storage.
//
// That is also the honest cost of a Doppler bank: coverage is bought in
// megabytes and in correlations per block. An HFM needing one bin where an LFM
// needs thirty is not a detail, it is the whole design argument.
template <std::size_t MaxTemplates, std::size_t FftSize>
class PulseBank {
    static_assert(MaxTemplates >= 1, "a bank needs at least one template");

  public:
    explicit PulseBank(Real sample_rate_hz, Real sound_speed_mps = 1500) noexcept
        : sample_rate_(sample_rate_hz), sound_speed_(sound_speed_mps) {}

    // Renders `spec` and prepares its matched filter. Returns false if the bank
    // is full, the spec is invalid, or the pulse is longer than the transform.
    bool add(const PulseSpec& spec) noexcept {
        if (count_ >= MaxTemplates) return false;
        if (!spec.valid()) return false;
        MatchedFilter mf;
        if (!matched_filter_from_pulse(fft_.view(), spec, sample_rate_,
                                       replica_scratch_,
                                       std::span<Complex>(spectra_[count_].data(), FftSize),
                                       mf)) {
            return false;
        }
        entries_[count_].spec = spec;
        entries_[count_].filter = mf;
        entries_[count_].doppler = 0;
        entries_[count_].doppler_bin = false;
        if (mf.replica_length > max_len_) max_len_ = mf.replica_length;
        ++count_;
        return true;
    }

    // Adds a replica time-scaled by (1 + doppler), i.e. matched to a
    // transmitter closing at v = doppler * c.
    bool add_at_doppler(const PulseSpec& spec, Real doppler, bool is_bin = true) noexcept {
        if (count_ >= MaxTemplates || !spec.valid()) return false;
        const std::size_t n = render_analytic_doppler(spec, sample_rate_, doppler,
                                                      replica_scratch_);
        if (n == 0) return false;
        MatchedFilter mf;
        if (!matched_filter_prepare(fft_.view(),
                                    std::span<const Complex>(replica_scratch_.data(), n),
                                    std::span<Complex>(spectra_[count_].data(), FftSize), mf)) {
            return false;
        }
        entries_[count_].spec = spec;
        entries_[count_].filter = mf;
        entries_[count_].doppler = doppler;
        entries_[count_].doppler_bin = is_bin;
        if (mf.replica_length > max_len_) max_len_ = mf.replica_length;
        ++count_;
        return true;
    }

    // Adds one waveform replicated across Doppler bins spanning
    // [v_min_mps, v_max_mps], spaced so the worst straddling loss between
    // adjacent bins stays under `max_loss_db`.
    //
    // How many bins that takes is a property of the waveform, not of the
    // library: an HFM of a given time-bandwidth product tolerates roughly two
    // orders of magnitude more Doppler than an LFM, so it may need a single
    // bin where the LFM needs dozens. Returns the number added, 0 on failure.
    // Nothing is added at all if the whole sweep would not fit.
    std::size_t add_doppler_bank(const PulseSpec& spec,
                                 Real v_min_mps, Real v_max_mps,
                                 Real max_loss_db = 1) noexcept {
        if (!spec.valid() || !(v_max_mps > v_min_mps)) return 0;
        // Straddling loss is worst midway between bins, so the spacing is twice
        // the tolerance at the requested loss.
        const Real tol = doppler_tolerance(spec, max_loss_db);
        if (!(tol > 0)) return 0;
        const Real step = 2 * tol * sound_speed_;

        // Round the interval count UP: with floor, the realised spacing can be
        // twice the design spacing, and the straddling-loss guarantee the
        // caller asked for silently does not hold.
        const Real span = v_max_mps - v_min_mps;
        auto intervals = static_cast<std::size_t>(span / step);
        if (static_cast<Real>(intervals) * step < span) ++intervals;
        const std::size_t bins = intervals + 1;
        if (count_ + bins > MaxTemplates) return 0;

        const Real actual_step = (bins > 1) ? span / static_cast<Real>(bins - 1) : 0;
        std::size_t added = 0;
        for (std::size_t i = 0; i < bins; ++i) {
            const Real v = v_min_mps + actual_step * static_cast<Real>(i);
            if (!add_at_doppler(spec, v / sound_speed_, true)) break;
            ++added;
        }
        return added;
    }

    // Bins the last add_doppler_bank() call would need, without adding them.
    [[nodiscard]] std::size_t doppler_bins_required(const PulseSpec& spec,
                                                    Real v_min_mps, Real v_max_mps,
                                                    Real max_loss_db = 1) const noexcept {
        if (!spec.valid() || !(v_max_mps > v_min_mps)) return 0;
        const Real tol = doppler_tolerance(spec, max_loss_db);
        if (!(tol > 0)) return 0;
        const Real step = 2 * tol * sound_speed_;
        const Real span = v_max_mps - v_min_mps;
        auto intervals = static_cast<std::size_t>(span / step);
        if (static_cast<Real>(intervals) * step < span) ++intervals;
        return intervals + 1;
    }

    // Velocity spacing the bins would actually be placed at, m/s. This is the
    // quantisation of any reported radial velocity, so a caller comparing an
    // estimate against truth should budget half of it.
    [[nodiscard]] Real doppler_bin_spacing_mps(const PulseSpec& spec,
                                               Real v_min_mps, Real v_max_mps,
                                               Real max_loss_db = 1) const noexcept {
        const std::size_t bins = doppler_bins_required(spec, v_min_mps, v_max_mps, max_loss_db);
        if (bins < 2) return 0;
        return (v_max_mps - v_min_mps) / static_cast<Real>(bins - 1);
    }

    void clear() noexcept { count_ = 0; max_len_ = 0; }

    [[nodiscard]] AnalyzerView view() const noexcept {
        return AnalyzerView{std::span<const PulseTemplate>(entries_.data(), count_),
                            fft_.view(), sample_rate_, max_len_, sound_speed_};
    }

    [[nodiscard]] std::size_t size() const noexcept { return count_; }
    [[nodiscard]] std::size_t max_replica_length() const noexcept { return max_len_; }
    [[nodiscard]] std::size_t stride() const noexcept {
        return matched_filter_stride(FftSize, max_len_);
    }
    [[nodiscard]] static constexpr std::size_t fft_size() noexcept { return FftSize; }

  private:
    FftPlan<FftSize> fft_;
    Real sample_rate_;
    Real sound_speed_;
    std::array<PulseTemplate, MaxTemplates> entries_{};
    std::array<std::array<Complex, FftSize>, MaxTemplates> spectra_{};
    std::array<Complex, FftSize> replica_scratch_{};
    std::size_t count_ = 0;
    std::size_t max_len_ = 0;
};

// Owning workspace to match a PulseBank of the same transform size.
template <std::size_t FftSize>
class AnalyzerScratch {
  public:
    [[nodiscard]] AnalyzerWorkspace view() noexcept {
        return AnalyzerWorkspace{
            std::span<Complex>(fft_scratch_.data(), FftSize),
            std::span<Complex>(correlation_.data(), FftSize),
            std::span<Real>(best_power_.data(), FftSize),
            std::span<std::uint16_t>(best_template_.data(), FftSize),
            std::span<Real>(prefix_.data(), FftSize + 1),
        };
    }

  private:
    std::array<Complex, FftSize>       fft_scratch_{};
    std::array<Complex, FftSize>       correlation_{};
    std::array<Real, FftSize>          best_power_{};
    std::array<std::uint16_t, FftSize> best_template_{};
    std::array<Real, FftSize + 1>      prefix_{};
};

}  // namespace phantom

#endif  // PHANTOM_PING_ANALYZER_HPP
