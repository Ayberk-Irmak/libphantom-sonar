/* SPDX-License-Identifier: Apache-2.0
 *
 * libphantom-sonar — C ABI.
 *
 * Hand-written, not generated. A generator would produce a faithful mirror of
 * the C++ headers, and a faithful mirror is the wrong thing: the C++ interface
 * uses spans, templates and RAII, none of which have a stable C representation,
 * and pretending otherwise is how an ABI breaks silently between releases.
 *
 * WHO OWNS THE MEMORY. The library allocates nothing, and this ABI does not
 * quietly change that by handing out pointers a C caller must free. Every
 * stateful object works the same way:
 *
 *     size_t n = ph_profile_size();
 *     void*  storage = malloc(n);        // or a static buffer, or a stack one
 *     ph_profile* p = ph_profile_init(storage, n);
 *
 * The caller decides where it lives, the library only decides how big it is and
 * what shape it has. That keeps a static-only or interrupt-context caller
 * possible, which is the whole reason the C++ side has the constraint.
 *
 * WHAT IS COVERED. Not everything. This exposes the subsystems that a caller
 * embedding the library actually reaches for -- sound speed, profiles and ray
 * tracing, waveform synthesis, tracking, and the communication codecs. The
 * beamformer, eigenray and reverberation paths are not here yet, because their
 * C++ interfaces take several spans at once and a C shape for them has not been
 * designed. "Stable C ABI" means what is here will not change shape, not that
 * everything is here.
 *
 * PRECISION. ph_real follows the build: double by default, float under
 * PHANTOM_REAL_FLOAT. A binary built one way cannot be linked against a caller
 * assuming the other, so ph_real_is_double() exists to check at runtime and
 * PH_REAL_IS_DOUBLE at compile time.
 *
 * THREADING. Nothing here holds global state and no function uses static
 * storage, so two threads may use two different objects freely. A single object
 * is not internally synchronised.
 *
 * STACK. ph_tracker_step converts its measurements through a stack buffer of
 * about 1.5 kB. Everything else is O(1) in stack. That matters on a part where
 * the whole stack is 8 kB.
 *
 * ERRORS. Functions that can fail return ph_status. Functions that cannot --
 * pure arithmetic on valid inputs -- return their value directly. There is no
 * errno-style global.
 */
#ifndef PHANTOM_PHANTOM_H
#define PHANTOM_PHANTOM_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ------------------------------------------------------------------------- */
/* Version and build configuration                                            */
/* ------------------------------------------------------------------------- */

#define PH_VERSION_MAJOR 0
#define PH_VERSION_MINOR 16
#define PH_VERSION_PATCH 0

/* Runtime version, in case a caller is linked against a different build than
 * it compiled against. Returns major*10000 + minor*100 + patch. */
int ph_version(void);
const char* ph_version_string(void);

#if defined(PHANTOM_REAL_FLOAT)
typedef float ph_real;
#define PH_REAL_IS_DOUBLE 0
#else
typedef double ph_real;
#define PH_REAL_IS_DOUBLE 1
#endif

/* Non-zero if the LIBRARY was built with double. Compare against
 * PH_REAL_IS_DOUBLE at startup: a mismatch means every ph_real crossing the
 * boundary is being reinterpreted, which produces plausible nonsense rather
 * than a crash. */
int ph_real_is_double(void);

typedef enum ph_status {
    PH_OK = 0,
    PH_ERR_NULL,          /* a required pointer was NULL */
    PH_ERR_SIZE,          /* a buffer was too small, or a size was zero */
    PH_ERR_RANGE,         /* an argument was outside its valid range */
    PH_ERR_STATE,         /* the object was not in a state allowing this call */
    PH_ERR_UNSUPPORTED    /* the build does not provide this */
} ph_status;

const char* ph_status_string(ph_status s);

/* ------------------------------------------------------------------------- */
/* Sound speed                                                                */
/* ------------------------------------------------------------------------- */

/* Depth in metres for the first three; Chen-Millero takes PRESSURE in bar.
 * See docs/math_spec.md §1 and §17 for which to use and why they differ. */
ph_real ph_sound_speed_medwin(ph_real temperature_c, ph_real salinity_psu, ph_real depth_m);
ph_real ph_sound_speed_mackenzie(ph_real temperature_c, ph_real salinity_psu, ph_real depth_m);

/* The original Chen & Millero (1977), IPTS-68. Reproduces the published UNESCO
 * check value; use it when comparing against pre-1990 tables. */
ph_real ph_sound_speed_chen_millero_1977(ph_real temperature_c, ph_real salinity_psu,
                                         ph_real pressure_bar);
/* Wong & Zhu (1995), ITS-90. The default, and what modern data wants. */
ph_real ph_sound_speed_chen_millero_its90(ph_real temperature_c, ph_real salinity_psu,
                                          ph_real pressure_bar);
/* Del Grosso (1974), ITS-90 form. Pressure in kg/cm^2. Independent of the
 * above -- different laboratory, different data, different functional form. */
ph_real ph_sound_speed_del_grosso(ph_real temperature_c, ph_real salinity_psu,
                                  ph_real pressure_kgcm2);

/* UNESCO from depth rather than pressure, via Leroy & Parthiot. */
ph_real ph_sound_speed_unesco(ph_real temperature_c, ph_real salinity_psu,
                              ph_real depth_m, ph_real latitude_deg);

ph_real ph_depth_to_pressure_bar(ph_real depth_m, ph_real latitude_deg);
ph_real ph_bar_to_kgcm2(ph_real bar);
ph_real ph_t90_to_t68(ph_real t90);

/* Munk's canonical deep sound channel, the standard ray-tracer benchmark. */
ph_real ph_sound_speed_munk(ph_real depth_m, ph_real axis_depth_m, ph_real axis_speed_mps,
                            ph_real epsilon, ph_real scale_m);

/* ------------------------------------------------------------------------- */
/* Sound speed profiles                                                       */
/* ------------------------------------------------------------------------- */

typedef struct ph_profile ph_profile;

/* Bytes of storage a profile needs, and the alignment it requires. */
size_t ph_profile_size(void);
size_t ph_profile_align(void);

/* Maximum number of (depth, speed) samples a profile holds. Fixed at build
 * time: this is a no-allocation library, so the capacity is a constant rather
 * than something the caller chooses. */
size_t ph_profile_capacity(void);

/* Places a profile in caller-supplied storage. `storage` must be at least
 * ph_profile_size() bytes and suitably aligned. Returns NULL on bad input.
 * There is no destructor: the object owns nothing. */
ph_profile* ph_profile_init(void* storage, size_t storage_bytes);

/* Appends a sample. Depths must strictly increase. */
ph_status ph_profile_push(ph_profile* p, ph_real depth_m, ph_real speed_mps);
void      ph_profile_clear(ph_profile* p);
size_t    ph_profile_count(const ph_profile* p);
ph_status ph_profile_sample(const ph_profile* p, size_t index,
                            ph_real* out_depth_m, ph_real* out_speed_mps);
/* Linear interpolation, clamped outside the profile. */
ph_real   ph_profile_speed_at(const ph_profile* p, ph_real depth_m);

/* ------------------------------------------------------------------------- */
/* Ray tracing                                                                */
/* ------------------------------------------------------------------------- */

typedef enum ph_boundary_action {
    PH_BOUNDARY_REFLECT = 0,
    PH_BOUNDARY_ABSORB  = 1
} ph_boundary_action;

typedef struct ph_trace_config {
    ph_real max_range_m;
    ph_real max_time_s;
    ph_real surface_depth_m;
    ph_real bottom_depth_m;      /* 0 means "the profile's deepest point" */
    ph_boundary_action surface;
    ph_boundary_action bottom;
    uint32_t max_bounces;
} ph_trace_config;

/* Fills `cfg` with the same defaults the C++ side uses. Call this rather than
 * zeroing the struct: a zeroed config traces zero range. */
void ph_trace_config_defaults(ph_trace_config* cfg);

typedef struct ph_ray_point {
    ph_real range_m;
    ph_real depth_m;
    ph_real angle_rad;     /* grazing angle, positive downgoing */
    ph_real time_s;
    ph_real speed_mps;
} ph_ray_point;

typedef enum ph_trace_status {
    PH_TRACE_MAX_RANGE = 0,
    PH_TRACE_MAX_TIME,
    PH_TRACE_BUFFER_FULL,
    PH_TRACE_ABSORBED,
    PH_TRACE_MAX_BOUNCES,
    PH_TRACE_DEGENERATE
} ph_trace_status;

typedef struct ph_trace_result {
    size_t   point_count;
    uint32_t surface_bounces;
    uint32_t bottom_bounces;
    uint32_t turning_points;
    uint32_t steps;
    ph_real  final_range_m;
    ph_real  final_time_s;
    ph_real  arc_length_m;
    ph_trace_status status;
} ph_trace_result;

/* Traces one ray. `out` receives up to `out_capacity` points.
 *
 * The trace writes directly into `out`; there is no internal buffer and no
 * bound on how many points a caller may ask for. Static assertions in the
 * implementation check that ph_ray_point and the C++ RayPoint agree in size,
 * alignment and every member offset, so this is verified at build time rather
 * than assumed. */
ph_status ph_trace_ray(const ph_profile* p,
                       ph_real source_depth_m,
                       ph_real launch_angle_rad,
                       const ph_trace_config* cfg,
                       ph_ray_point* out, size_t out_capacity,
                       ph_trace_result* result);

/* ------------------------------------------------------------------------- */
/* Waveforms                                                                  */
/* ------------------------------------------------------------------------- */

typedef enum ph_pulse_type {
    PH_PULSE_UNKNOWN = 0,
    PH_PULSE_CW,
    PH_PULSE_LFM_UP,
    PH_PULSE_LFM_DOWN,
    PH_PULSE_HFM
} ph_pulse_type;

typedef enum ph_taper {
    PH_TAPER_RECTANGULAR = 0,
    PH_TAPER_HANN,
    PH_TAPER_TUKEY25
} ph_taper;

typedef struct ph_pulse_spec {
    ph_pulse_type type;
    ph_real f_start_hz;
    ph_real f_end_hz;
    ph_real duration_s;
    ph_real amplitude;
    ph_taper taper;
} ph_pulse_spec;

/* Samples the pulse would occupy at this sample rate. */
size_t ph_pulse_length(const ph_pulse_spec* spec, ph_real sample_rate_hz);

/* Renders the real waveform. Returns samples written via `out_written`. */
ph_status ph_render_real(const ph_pulse_spec* spec, ph_real sample_rate_hz,
                         ph_real* out, size_t out_capacity, size_t* out_written);

/* Renders with the time axis scaled by (1 + doppler), i.e. what a receiver sees
 * from a target closing at v = doppler * c. NOTE the argument is v/c, not
 * 1 + v/c -- getting that wrong is a mistake this project has made. */
ph_status ph_render_real_doppler(const ph_pulse_spec* spec, ph_real sample_rate_hz,
                                 ph_real doppler,
                                 ph_real* out, size_t out_capacity, size_t* out_written);

/* ------------------------------------------------------------------------- */
/* Tracking                                                                   */
/* ------------------------------------------------------------------------- */

typedef struct ph_measurement {
    ph_real range_m;
    ph_real bearing_rad;
    ph_real time_s;
    ph_real range_rate_mps;
    int     has_range_rate;      /* 0 or 1; see the C++ header for why this is
                                  * separate from range_rate_mps being zero */
    uint16_t label;
    ph_real amplitude;
} ph_measurement;

typedef struct ph_tracker_config {
    ph_real process_accel_mps2;
    ph_real range_sigma_m;
    ph_real bearing_sigma_rad;
    ph_real range_rate_sigma_mps;
    ph_real gate_chi2_2dof;
    ph_real gate_chi2_3dof;
    uint32_t confirm_hits;
    uint32_t delete_misses;
    ph_real init_velocity_sigma_mps;
} ph_tracker_config;

void ph_tracker_config_defaults(ph_tracker_config* cfg);

typedef enum ph_track_status {
    PH_TRACK_FREE = 0,
    PH_TRACK_TENTATIVE,
    PH_TRACK_CONFIRMED,
    PH_TRACK_COASTING
} ph_track_status;

typedef struct ph_track_state {
    ph_real x, y, vx, vy;        /* y is along broadside; bearing = atan2(x, y) */
    ph_real range_m;
    ph_real bearing_rad;
    ph_real range_rate_mps;
    uint32_t id;
    uint32_t hits;
    uint32_t misses;
    ph_track_status status;
} ph_track_state;

typedef struct ph_tracker ph_tracker;

/* Storage for a tracker holding `max_tracks` tracks. */
size_t ph_tracker_size(size_t max_tracks);
size_t ph_tracker_align(void);
size_t ph_tracker_max_tracks(void);   /* the largest max_tracks accepted */

ph_tracker* ph_tracker_init(void* storage, size_t storage_bytes, size_t max_tracks);

/* One scan. Returns the number of live tracks via `out_live`. */
ph_status ph_tracker_step(ph_tracker* t,
                          const ph_measurement* measurements, size_t n_measurements,
                          const ph_tracker_config* cfg, ph_real time_s,
                          size_t* out_live);

/* Reads back track `index` in [0, max_tracks). */
ph_status ph_tracker_get(const ph_tracker* t, size_t index, ph_track_state* out);

/* Confirmed plus coasting -- the tracks a caller would call real. */
size_t ph_tracker_established(const ph_tracker* t);

/* Chi-square gate for a probability and 1, 2 or 3 degrees of freedom. */
ph_real ph_chi2_gate(ph_real probability, size_t dof);

/* ------------------------------------------------------------------------- */
/* Communication                                                              */
/* ------------------------------------------------------------------------- */

/* Maximal-length PN sequence, 2^degree - 1 chips of +1/-1. Degree 5..15. */
size_t ph_msequence_length(uint32_t degree);
ph_status ph_generate_msequence(uint32_t degree, uint32_t seed,
                                ph_real* out, size_t out_capacity, size_t* out_written);

/* 10*log10(N). Read the C++ header before using this number: against white
 * noise at a fixed energy per bit, spreading buys nothing. */
ph_real ph_processing_gain_db(size_t chips_per_bit);

/* Chip slip across one bit at a given closing speed. Once this approaches half
 * a chip the coherent sum has cancelled itself. */
ph_real ph_chip_slip(size_t chips_per_bit, ph_real closing_speed_mps,
                     ph_real sound_speed_mps);

uint32_t ph_crc32(const uint8_t* data, size_t length);

/* Reed-Solomon (15, 11) over GF(16): 11 data symbols in, 15 out, each < 16.
 * Corrects up to 2 symbol errors. */
#define PH_RS_N 15
#define PH_RS_K 11

ph_status ph_rs_encode(const uint8_t* data11, uint8_t* out15);

typedef enum ph_rs_result {
    PH_RS_CLEAN = 0,
    PH_RS_CORRECTED,
    PH_RS_UNCORRECTABLE,
    PH_RS_BAD_INPUT
} ph_rs_result;

/* Decodes in place. `out_corrected` may be NULL. */
ph_rs_result ph_rs_decode(uint8_t* codeword15, size_t* out_corrected);

#ifdef __cplusplus
}  /* extern "C" */
#endif

#endif /* PHANTOM_PHANTOM_H */
