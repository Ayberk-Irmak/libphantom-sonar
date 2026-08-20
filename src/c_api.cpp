// SPDX-License-Identifier: Apache-2.0
//
// The C ABI, implemented over the C++ interface.
//
// Two rules hold throughout and are worth stating once rather than repeating at
// every function:
//
//   NOTHING IS ALLOCATED. Every stateful type is placement-constructed into
//   storage the caller supplied. The `_size()` and `_align()` queries exist so
//   the caller can provide it from wherever suits -- static, stack, or a pool --
//   and the library never learns which.
//
//   NOTHING THROWS. The C++ side is noexcept end to end, so there is no
//   try/catch here. If that ever stops being true, this file is where the
//   breakage would escape into C, and the -fno-exceptions build mode added in
//   this release is what stops it silently.
#include "phantom/phantom.h"

#include "phantom/comm.hpp"
#include "phantom/config.hpp"
#include "phantom/profile.hpp"
#include "phantom/ray_tracer.hpp"
#include "phantom/sound_speed.hpp"
#include "phantom/tracker.hpp"
#include "phantom/types.hpp"
#include "phantom/waveform.hpp"

#include <cstddef>   // offsetof, for the ray-point layout assertions
#include <new>       // placement new only; nothing here allocates
#include <type_traits>

namespace {

// The one capacity the C ABI freezes. Changing it is an ABI break, which is why
// it is named here rather than left implicit in a template argument scattered
// through the file.
constexpr std::size_t kProfileCapacity = 2048;
constexpr std::size_t kMaxTracks = 64;

using CProfile = phantom::SoundSpeedProfile<kProfileCapacity>;

struct TrackerImpl {
    // Default member initialisers, like every other struct in the library.
    // ph_tracker_init placement-constructs this and then sets both fields, so
    // they were never read uninitialised -- but relying on that is a promise
    // about a call site rather than a property of the type.
    std::uint32_t next_id = 1;
    std::size_t   count = 0;
    phantom::Track tracks[kMaxTracks] = {};
};

phantom::Real R(ph_real v) noexcept { return static_cast<phantom::Real>(v); }
ph_real C(phantom::Real v) noexcept { return static_cast<ph_real>(v); }

phantom::TraceConfig to_cpp(const ph_trace_config& c) noexcept {
    phantom::TraceConfig out;
    out.max_range_m     = R(c.max_range_m);
    out.max_time_s      = R(c.max_time_s);
    out.surface_depth_m = R(c.surface_depth_m);
    out.bottom_depth_m  = R(c.bottom_depth_m);
    out.surface = (c.surface == PH_BOUNDARY_ABSORB) ? phantom::BoundaryAction::Absorb
                                                    : phantom::BoundaryAction::Reflect;
    out.bottom  = (c.bottom == PH_BOUNDARY_ABSORB) ? phantom::BoundaryAction::Absorb
                                                   : phantom::BoundaryAction::Reflect;
    out.max_bounces = c.max_bounces;
    return out;
}

phantom::PulseSpec to_cpp(const ph_pulse_spec& s) noexcept {
    phantom::PulseSpec out;
    switch (s.type) {
        case PH_PULSE_CW:       out.type = phantom::PulseType::Cw; break;
        case PH_PULSE_LFM_UP:   out.type = phantom::PulseType::LfmUp; break;
        case PH_PULSE_LFM_DOWN: out.type = phantom::PulseType::LfmDown; break;
        case PH_PULSE_HFM:      out.type = phantom::PulseType::Hfm; break;
        default:                out.type = phantom::PulseType::Unknown; break;
    }
    out.f_start_hz = R(s.f_start_hz);
    out.f_end_hz   = R(s.f_end_hz);
    out.duration_s = R(s.duration_s);
    out.amplitude  = R(s.amplitude);
    switch (s.taper) {
        case PH_TAPER_HANN:    out.taper = phantom::Taper::Hann; break;
        case PH_TAPER_TUKEY25: out.taper = phantom::Taper::Tukey25; break;
        default:               out.taper = phantom::Taper::Rectangular; break;
    }
    return out;
}

phantom::TrackerConfig to_cpp(const ph_tracker_config& c) noexcept {
    phantom::TrackerConfig out;
    out.process_accel_mps2      = R(c.process_accel_mps2);
    out.range_sigma_m           = R(c.range_sigma_m);
    out.bearing_sigma_rad       = R(c.bearing_sigma_rad);
    out.range_rate_sigma_mps    = R(c.range_rate_sigma_mps);
    out.gate_chi2_2dof          = R(c.gate_chi2_2dof);
    out.gate_chi2_3dof          = R(c.gate_chi2_3dof);
    out.confirm_hits            = c.confirm_hits;
    out.delete_misses           = c.delete_misses;
    out.init_velocity_sigma_mps = R(c.init_velocity_sigma_mps);
    return out;
}

phantom::Measurement to_cpp(const ph_measurement& m) noexcept {
    phantom::Measurement out;
    out.range_m        = R(m.range_m);
    out.bearing_rad    = R(m.bearing_rad);
    out.time_s         = R(m.time_s);
    out.range_rate_mps = R(m.range_rate_mps);
    out.has_range_rate = (m.has_range_rate != 0);
    out.label          = m.label;
    out.amplitude      = R(m.amplitude);
    return out;
}

}  // namespace

extern "C" {

// --- version -----------------------------------------------------------------

int ph_version(void) {
    return PHANTOM_VERSION_MAJOR * 10000 + PHANTOM_VERSION_MINOR * 100 + PHANTOM_VERSION_PATCH;
}

const char* ph_version_string(void) { return PHANTOM_VERSION_STRING; }

int ph_real_is_double(void) {
    // Not sizeof(Real) == sizeof(double): in a double build that is literally
    // sizeof(double) == sizeof(double), which constant-folds to a tautology and
    // says nothing about the type. This asks the question that is meant.
    return std::is_same_v<phantom::Real, double> ? 1 : 0;
}

const char* ph_status_string(ph_status s) {
    switch (s) {
        case PH_OK:              return "ok";
        case PH_ERR_NULL:        return "null pointer";
        case PH_ERR_SIZE:        return "buffer too small or size zero";
        case PH_ERR_RANGE:       return "argument out of range";
        case PH_ERR_STATE:       return "object not in a state allowing this call";
        case PH_ERR_UNSUPPORTED: return "not supported by this build";
    }
    return "unknown";
}

// --- sound speed -------------------------------------------------------------

ph_real ph_sound_speed_medwin(ph_real t, ph_real s, ph_real z) {
    return C(phantom::sound_speed::medwin(R(t), R(s), R(z)));
}
ph_real ph_sound_speed_mackenzie(ph_real t, ph_real s, ph_real z) {
    return C(phantom::sound_speed::mackenzie(R(t), R(s), R(z)));
}
ph_real ph_sound_speed_chen_millero_1977(ph_real t, ph_real s, ph_real p) {
    return C(phantom::sound_speed::chen_millero_1977(R(t), R(s), R(p)));
}
ph_real ph_sound_speed_chen_millero_its90(ph_real t, ph_real s, ph_real p) {
    return C(phantom::sound_speed::chen_millero_its90(R(t), R(s), R(p)));
}
ph_real ph_sound_speed_del_grosso(ph_real t, ph_real s, ph_real p) {
    return C(phantom::sound_speed::del_grosso(R(t), R(s), R(p)));
}
ph_real ph_sound_speed_unesco(ph_real t, ph_real s, ph_real z, ph_real lat) {
    return C(phantom::sound_speed::unesco(R(t), R(s), R(z), R(lat)));
}
ph_real ph_depth_to_pressure_bar(ph_real z, ph_real lat) {
    return C(phantom::sound_speed::depth_to_pressure_bar(R(z), R(lat)));
}
ph_real ph_bar_to_kgcm2(ph_real bar) { return C(phantom::sound_speed::bar_to_kgcm2(R(bar))); }
ph_real ph_t90_to_t68(ph_real t90) { return C(phantom::sound_speed::t90_to_t68(R(t90))); }

ph_real ph_sound_speed_munk(ph_real z, ph_real axis_z, ph_real axis_c,
                            ph_real eps, ph_real scale) {
    return C(phantom::sound_speed::munk(R(z), R(axis_z), R(axis_c), R(eps), R(scale)));
}

// --- profiles ----------------------------------------------------------------

size_t ph_profile_size(void)     { return sizeof(CProfile); }
size_t ph_profile_align(void)    { return alignof(CProfile); }
size_t ph_profile_capacity(void) { return kProfileCapacity; }

ph_profile* ph_profile_init(void* storage, size_t storage_bytes) {
    if (storage == nullptr || storage_bytes < sizeof(CProfile)) return nullptr;
    const auto addr = reinterpret_cast<std::uintptr_t>(storage);
    if (addr % alignof(CProfile) != 0) return nullptr;
    return reinterpret_cast<ph_profile*>(new (storage) CProfile());
}

ph_status ph_profile_push(ph_profile* p, ph_real depth_m, ph_real speed_mps) {
    if (p == nullptr) return PH_ERR_NULL;
    auto* impl = reinterpret_cast<CProfile*>(p);
    return impl->push(R(depth_m), R(speed_mps)) ? PH_OK : PH_ERR_RANGE;
}

void ph_profile_clear(ph_profile* p) {
    if (p != nullptr) reinterpret_cast<CProfile*>(p)->clear();
}

size_t ph_profile_count(const ph_profile* p) {
    return (p == nullptr) ? 0 : reinterpret_cast<const CProfile*>(p)->size();
}

ph_status ph_profile_sample(const ph_profile* p, size_t index,
                            ph_real* out_depth_m, ph_real* out_speed_mps) {
    if (p == nullptr) return PH_ERR_NULL;
    const auto* impl = reinterpret_cast<const CProfile*>(p);
    if (index >= impl->size()) return PH_ERR_RANGE;
    const phantom::ProfileView v = impl->view();
    if (out_depth_m != nullptr) *out_depth_m = C(v.depth_m[index]);
    if (out_speed_mps != nullptr) *out_speed_mps = C(v.speed_mps[index]);
    return PH_OK;
}

ph_real ph_profile_speed_at(const ph_profile* p, ph_real depth_m) {
    if (p == nullptr) return 0;
    const auto* impl = reinterpret_cast<const CProfile*>(p);
    if (!impl->valid()) return 0;
    return C(phantom::speed_at(impl->view(), R(depth_m)));
}

// --- ray tracing -------------------------------------------------------------

void ph_trace_config_defaults(ph_trace_config* cfg) {
    if (cfg == nullptr) return;
    const phantom::TraceConfig d;
    cfg->max_range_m     = C(d.max_range_m);
    cfg->max_time_s      = C(d.max_time_s);
    cfg->surface_depth_m = C(d.surface_depth_m);
    cfg->bottom_depth_m  = C(d.bottom_depth_m);
    cfg->surface = PH_BOUNDARY_REFLECT;
    cfg->bottom  = PH_BOUNDARY_REFLECT;
    cfg->max_bounces = d.max_bounces;
}

ph_status ph_trace_ray(const ph_profile* p, ph_real source_depth_m, ph_real launch_angle_rad,
                       const ph_trace_config* cfg,
                       ph_ray_point* out, size_t out_capacity, ph_trace_result* result) {
    if (p == nullptr || cfg == nullptr || out == nullptr || result == nullptr) return PH_ERR_NULL;
    if (out_capacity == 0) return PH_ERR_SIZE;
    const auto* impl = reinterpret_cast<const CProfile*>(p);
    if (!impl->valid()) return PH_ERR_STATE;

    // Trace DIRECTLY into the caller's array.
    //
    // The first version copied through a static 8192-point scratch buffer,
    // "because the ABI must not assume the two structs share a layout". That
    // was 160 kB of .bss -- half the RAM of the Cortex-M7 this library is meant
    // to fit on -- and it made the function non-reentrant, both for a copy that
    // does nothing. The cross-compilation size report is what exposed it.
    //
    // The assumption can be CHECKED instead of avoided. If any of these ever
    // stops holding, the build breaks here rather than the ABI silently
    // scrambling ray paths.
    static_assert(sizeof(ph_ray_point) == sizeof(phantom::RayPoint),
                  "ph_ray_point and RayPoint must agree in size");
    static_assert(alignof(ph_ray_point) == alignof(phantom::RayPoint),
                  "ph_ray_point and RayPoint must agree in alignment");
    static_assert(offsetof(ph_ray_point, range_m)   == offsetof(phantom::RayPoint, range_m), "");
    static_assert(offsetof(ph_ray_point, depth_m)   == offsetof(phantom::RayPoint, depth_m), "");
    static_assert(offsetof(ph_ray_point, angle_rad) == offsetof(phantom::RayPoint, angle_rad), "");
    static_assert(offsetof(ph_ray_point, time_s)    == offsetof(phantom::RayPoint, time_s), "");
    static_assert(offsetof(ph_ray_point, speed_mps) == offsetof(phantom::RayPoint, speed_mps), "");
    static_assert(sizeof(ph_real) == sizeof(phantom::Real),
                  "ph_real must match Real; check PHANTOM_REAL_FLOAT on both sides");

    const phantom::TraceConfig c = to_cpp(*cfg);
    const phantom::TraceResult r =
        phantom::trace_ray(impl->view(), R(source_depth_m), R(launch_angle_rad), c,
                           std::span<phantom::RayPoint>(
                               reinterpret_cast<phantom::RayPoint*>(out), out_capacity));

    result->point_count     = r.point_count;
    result->surface_bounces = r.surface_bounces;
    result->bottom_bounces  = r.bottom_bounces;
    result->turning_points  = r.turning_points;
    result->steps           = r.steps;
    result->final_range_m   = C(r.final_range_m);
    result->final_time_s    = C(r.final_time_s);
    result->arc_length_m    = C(r.path_length_m);
    result->status          = static_cast<ph_trace_status>(static_cast<int>(r.status));
    return PH_OK;
}

// --- waveforms ---------------------------------------------------------------

size_t ph_pulse_length(const ph_pulse_spec* spec, ph_real fs) {
    if (spec == nullptr) return 0;
    const phantom::PulseSpec s = to_cpp(*spec);
    return phantom::pulse_length(s, R(fs));
}

ph_status ph_render_real(const ph_pulse_spec* spec, ph_real fs,
                         ph_real* out, size_t out_capacity, size_t* out_written) {
    if (spec == nullptr || out == nullptr) return PH_ERR_NULL;
    const phantom::PulseSpec s = to_cpp(*spec);
    // ph_real and phantom::Real are the same type by construction; the C header
    // derives one from the same PHANTOM_REAL_FLOAT that picks the other, and
    // ph_trace_ray static_asserts it. An earlier version guarded this with an
    // #if whose two branches were identical, which guarded nothing.
    const std::size_t n = phantom::render_real(s, R(fs),
        std::span<phantom::Real>(reinterpret_cast<phantom::Real*>(out), out_capacity));
    if (out_written != nullptr) *out_written = n;
    return (n > 0) ? PH_OK : PH_ERR_SIZE;
}

ph_status ph_render_real_doppler(const ph_pulse_spec* spec, ph_real fs, ph_real doppler,
                                 ph_real* out, size_t out_capacity, size_t* out_written) {
    if (spec == nullptr || out == nullptr) return PH_ERR_NULL;
    const phantom::PulseSpec s = to_cpp(*spec);
    const std::size_t n = phantom::render_real_doppler(s, R(fs), R(doppler),
        std::span<phantom::Real>(reinterpret_cast<phantom::Real*>(out), out_capacity));
    if (out_written != nullptr) *out_written = n;
    return (n > 0) ? PH_OK : PH_ERR_SIZE;
}

// --- tracking ----------------------------------------------------------------

void ph_tracker_config_defaults(ph_tracker_config* cfg) {
    if (cfg == nullptr) return;
    const phantom::TrackerConfig d;
    cfg->process_accel_mps2      = C(d.process_accel_mps2);
    cfg->range_sigma_m           = C(d.range_sigma_m);
    cfg->bearing_sigma_rad       = C(d.bearing_sigma_rad);
    cfg->range_rate_sigma_mps    = C(d.range_rate_sigma_mps);
    cfg->gate_chi2_2dof          = C(d.gate_chi2_2dof);
    cfg->gate_chi2_3dof          = C(d.gate_chi2_3dof);
    cfg->confirm_hits            = d.confirm_hits;
    cfg->delete_misses           = d.delete_misses;
    cfg->init_velocity_sigma_mps = C(d.init_velocity_sigma_mps);
}

size_t ph_tracker_size(size_t max_tracks) {
    if (max_tracks == 0 || max_tracks > kMaxTracks) return 0;
    return sizeof(TrackerImpl);
}
size_t ph_tracker_align(void) { return alignof(TrackerImpl); }
size_t ph_tracker_max_tracks(void) { return kMaxTracks; }

ph_tracker* ph_tracker_init(void* storage, size_t storage_bytes, size_t max_tracks) {
    if (storage == nullptr || max_tracks == 0 || max_tracks > kMaxTracks) return nullptr;
    if (storage_bytes < sizeof(TrackerImpl)) return nullptr;
    if (reinterpret_cast<std::uintptr_t>(storage) % alignof(TrackerImpl) != 0) return nullptr;
    auto* impl = new (storage) TrackerImpl();
    impl->next_id = 1;
    impl->count = max_tracks;
    for (std::size_t i = 0; i < kMaxTracks; ++i) impl->tracks[i] = phantom::Track{};
    return reinterpret_cast<ph_tracker*>(impl);
}

ph_status ph_tracker_step(ph_tracker* t, const ph_measurement* measurements,
                          size_t n_measurements, const ph_tracker_config* cfg,
                          ph_real time_s, size_t* out_live) {
    if (t == nullptr || cfg == nullptr) return PH_ERR_NULL;
    if (n_measurements > 0 && measurements == nullptr) return PH_ERR_NULL;
    auto* impl = reinterpret_cast<TrackerImpl*>(t);

    // A stack buffer, not a static one. ph_measurement and phantom::Measurement
    // genuinely differ -- the C side carries an int where C++ has a bool -- so a
    // conversion is unavoidable here, unlike the ray-point path above. Putting
    // it on the stack costs about 1.5 kB of stack and buys reentrancy, which is
    // the right way round for a library that may be called from two threads.
    phantom::Measurement scratch[64];
    const std::size_t n = (n_measurements < 64) ? n_measurements : std::size_t{64};
    for (std::size_t i = 0; i < n; ++i) scratch[i] = to_cpp(measurements[i]);

    const phantom::TrackerConfig c = to_cpp(*cfg);
    const std::size_t live = phantom::tracker_step(
        std::span<phantom::Track>(impl->tracks, impl->count),
        std::span<const phantom::Measurement>(scratch, n), c, R(time_s), impl->next_id);
    if (out_live != nullptr) *out_live = live;
    return PH_OK;
}

ph_status ph_tracker_get(const ph_tracker* t, size_t index, ph_track_state* out) {
    if (t == nullptr || out == nullptr) return PH_ERR_NULL;
    const auto* impl = reinterpret_cast<const TrackerImpl*>(t);
    if (index >= impl->count) return PH_ERR_RANGE;
    const phantom::Track& tr = impl->tracks[index];
    out->x  = C(tr.state.x);
    out->y  = C(tr.state.y);
    out->vx = C(tr.state.vx);
    out->vy = C(tr.state.vy);
    out->range_m        = C(tr.state.range_m());
    out->bearing_rad    = C(tr.state.bearing_rad());
    out->range_rate_mps = C(tr.state.range_rate_mps());
    out->id     = tr.id;
    out->hits   = tr.hits;
    out->misses = tr.misses;
    out->status = static_cast<ph_track_status>(static_cast<int>(tr.status));
    return PH_OK;
}

size_t ph_tracker_established(const ph_tracker* t) {
    if (t == nullptr) return 0;
    const auto* impl = reinterpret_cast<const TrackerImpl*>(t);
    return phantom::count_established(
        std::span<const phantom::Track>(impl->tracks, impl->count));
}

ph_real ph_chi2_gate(ph_real probability, size_t dof) {
    return C(phantom::chi2_gate(R(probability), dof));
}

// --- communication -----------------------------------------------------------

size_t ph_msequence_length(uint32_t degree) {
    return phantom::comm::msequence_length(degree);
}

ph_status ph_generate_msequence(uint32_t degree, uint32_t seed,
                                ph_real* out, size_t out_capacity, size_t* out_written) {
    if (out == nullptr) return PH_ERR_NULL;
    const std::size_t n = phantom::comm::generate_msequence(degree, seed,
        std::span<phantom::Real>(reinterpret_cast<phantom::Real*>(out), out_capacity));
    if (out_written != nullptr) *out_written = n;
    if (n == 0) {
        return (phantom::comm::msequence_length(degree) == 0) ? PH_ERR_RANGE : PH_ERR_SIZE;
    }
    return PH_OK;
}

ph_real ph_processing_gain_db(size_t chips_per_bit) {
    return C(phantom::comm::processing_gain_db(chips_per_bit));
}

ph_real ph_chip_slip(size_t chips_per_bit, ph_real closing_speed_mps, ph_real c_mps) {
    phantom::comm::DsssConfig cfg;
    cfg.chips_per_bit = chips_per_bit;
    return C(phantom::comm::chip_slip(cfg, R(closing_speed_mps), R(c_mps)));
}

uint32_t ph_crc32(const uint8_t* data, size_t length) {
    if (data == nullptr && length > 0) return 0;
    return phantom::comm::crc32(std::span<const std::uint8_t>(data, length));
}

ph_status ph_rs_encode(const uint8_t* data11, uint8_t* out15) {
    if (data11 == nullptr || out15 == nullptr) return PH_ERR_NULL;
    const bool ok = phantom::comm::rs_encode(
        std::span<const std::uint8_t>(data11, PH_RS_K),
        std::span<std::uint8_t>(out15, PH_RS_N));
    return ok ? PH_OK : PH_ERR_RANGE;
}

ph_rs_result ph_rs_decode(uint8_t* codeword15, size_t* out_corrected) {
    if (codeword15 == nullptr) return PH_RS_BAD_INPUT;
    std::size_t fixed = 0;
    const phantom::comm::RsResult r = phantom::comm::rs_decode(
        std::span<std::uint8_t>(codeword15, PH_RS_N), fixed);
    if (out_corrected != nullptr) *out_corrected = fixed;
    return static_cast<ph_rs_result>(static_cast<int>(r));
}

}  // extern "C"
