// SPDX-License-Identifier: Apache-2.0
#include "phantom/reverberation.hpp"

#include <cmath>

namespace phantom {
namespace {

constexpr Real kZero = static_cast<Real>(0);
constexpr Real kTwo  = static_cast<Real>(2);

// A scattering strength below this is indistinguishable from "nothing scatters"
// and the logarithm of a zero sine would otherwise run to negative infinity.
constexpr Real kMinScatteringDb = static_cast<Real>(-120);

}  // namespace

Real lambert_bottom_scattering_db(Real mu_db, Real grazing_angle_rad) noexcept {
    const Real th = std::fabs(clamp(grazing_angle_rad, -kHalfPi, kHalfPi));
    const Real s = std::sin(th);
    if (!(s > kZero)) return kMinScatteringDb;
    const Real v = mu_db + static_cast<Real>(20) * std::log10(s);
    return (v < kMinScatteringDb) ? kMinScatteringDb : v;
}

Real chapman_harris_surface_scattering_db(Real wind_speed_knots,
                                          Real frequency_hz,
                                          Real grazing_angle_rad) noexcept {
    if (!(wind_speed_knots > kZero) || !(frequency_hz > kZero)) return kMinScatteringDb;
    const Real th_deg = rad2deg(std::fabs(clamp(grazing_angle_rad, -kHalfPi, kHalfPi)));
    if (!(th_deg > kZero)) return kMinScatteringDb;

    const Real beta = static_cast<Real>(158)
                    * std::pow(wind_speed_knots * std::cbrt(frequency_hz),
                               static_cast<Real>(-0.58));
    if (!(beta > kZero)) return kMinScatteringDb;

    const Real v = static_cast<Real>(3.3) * beta * std::log10(th_deg / static_cast<Real>(30))
                 - static_cast<Real>(42.4) * std::log10(beta)
                 + static_cast<Real>(2.6);
    return (v < kMinScatteringDb) ? kMinScatteringDb : v;
}

Real ensonified_area_m2(Real range_m, Real beamwidth_rad,
                        Real pulse_length_s, Real sound_speed_mps) noexcept {
    if (!(range_m > kZero) || !(beamwidth_rad > kZero)) return kZero;
    if (!(pulse_length_s > kZero) || !(sound_speed_mps > kZero)) return kZero;
    return range_m * beamwidth_rad * (sound_speed_mps * pulse_length_s / kTwo);
}

Real ensonified_volume_m3(Real range_m, Real solid_angle_sr,
                          Real pulse_length_s, Real sound_speed_mps) noexcept {
    if (!(range_m > kZero) || !(solid_angle_sr > kZero)) return kZero;
    if (!(pulse_length_s > kZero) || !(sound_speed_mps > kZero)) return kZero;
    return range_m * range_m * solid_angle_sr * (sound_speed_mps * pulse_length_s / kTwo);
}

Real reverberation_level_db(Real source_level_db,
                            Real transmission_loss_db,
                            Real scattering_strength_db,
                            Real scatterer_extent) noexcept {
    if (!(scatterer_extent > kZero)) return kMinScatteringDb;
    return source_level_db - kTwo * transmission_loss_db + scattering_strength_db
         + static_cast<Real>(10) * std::log10(scatterer_extent);
}

Real echo_to_reverberation_ratio_db(Real target_strength_db,
                                    Real scattering_strength_db,
                                    Real scatterer_extent) noexcept {
    if (!(scatterer_extent > kZero)) return kZero;
    // Source level and transmission loss are absent because they cancel; see
    // the derivation at the top of the header.
    return target_strength_db - scattering_strength_db
         - static_cast<Real>(10) * std::log10(scatterer_extent);
}

Real reverberation_limited_range_m(Real source_level_db,
                                   Real scattering_strength_db,
                                   Real beamwidth_rad,
                                   Real pulse_length_s,
                                   Real sound_speed_mps,
                                   Real noise_level_db,
                                   Real min_range_m,
                                   Real max_range_m) noexcept {
    if (!(max_range_m > min_range_m) || !(min_range_m > kZero)) return kZero;

    auto excess = [&](Real r) noexcept {
        const Real area = ensonified_area_m2(r, beamwidth_rad, pulse_length_s, sound_speed_mps);
        const Real tl = static_cast<Real>(20) * std::log10(r);   // spherical
        return reverberation_level_db(source_level_db, tl, scattering_strength_db, area)
             - noise_level_db;
    };

    // Reverberation falls monotonically as 30 log10(r), so if it is already
    // below the noise at the near edge the geometry is never reverb-limited.
    if (!(excess(min_range_m) > kZero)) return kZero;
    if (excess(max_range_m) > kZero) return max_range_m;

    Real lo = min_range_m;
    Real hi = max_range_m;
    for (int i = 0; i < 80; ++i) {
        const Real mid = (lo + hi) / kTwo;
        if (excess(mid) > kZero) lo = mid; else hi = mid;
    }
    return lo;
}

std::size_t reverberation_envelope(const ReverbProfile& profile,
                                   Real sample_rate_hz,
                                   Real start_time_s,
                                   std::span<Real> out) noexcept {
    if (!(sample_rate_hz > kZero) || out.empty()) return 0;
    if (!(profile.sound_speed_mps > kZero)) return 0;

    const Real dt = static_cast<Real>(1) / sample_rate_hz;
    const Real noise_amp = std::pow(static_cast<Real>(10),
                                    (profile.noise_level_db - profile.reference_level_db)
                                    / static_cast<Real>(20));

    for (std::size_t i = 0; i < out.size(); ++i) {
        const Real t = start_time_s + static_cast<Real>(i) * dt;
        const Real r = profile.sound_speed_mps * t / kTwo;
        if (!(r > kZero)) {
            out[i] = noise_amp;
            continue;
        }
        const Real area = ensonified_area_m2(r, profile.beamwidth_rad,
                                             profile.pulse_length_s, profile.sound_speed_mps);
        const Real tl = static_cast<Real>(20) * std::log10(r);
        const Real rl = reverberation_level_db(profile.source_level_db, tl,
                                               profile.scattering_strength_db, area);
        const Real amp = std::pow(static_cast<Real>(10),
                                  (rl - profile.reference_level_db) / static_cast<Real>(20));
        // The ambient floor does not decay, so it takes over once the
        // reverberation has fallen through it. Powers add, not amplitudes.
        out[i] = std::sqrt(amp * amp + noise_amp * noise_amp);
    }
    return out.size();
}

}  // namespace phantom
