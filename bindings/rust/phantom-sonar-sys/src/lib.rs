// SPDX-License-Identifier: Apache-2.0
//! Raw FFI declarations for `libphantom-sonar`, transcribed by hand from
//! `include/phantom/phantom.h`.
//!
//! Everything here is `unsafe` to call and carries no invariants. The safe
//! wrapper lives in the `phantom-sonar` crate; this one exists so that crate
//! has something to be safe *about*.
//!
//! # Precision
//!
//! `PhReal` is `f64` unless the C library was built with `PHANTOM_REAL_FLOAT`,
//! in which case it is `f32` and this crate must be built with the
//! `real_float` feature to match. A mismatch is not a link error -- it silently
//! reinterprets every number crossing the boundary -- so
//! [`ph_real_is_double`] exists and the safe wrapper checks it at construction.
#![allow(non_camel_case_types)]

use core::ffi::{c_char, c_int};

#[cfg(feature = "real_float")]
pub type PhReal = f32;
#[cfg(not(feature = "real_float"))]
pub type PhReal = f64;

#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PhStatus {
    Ok = 0,
    ErrNull = 1,
    ErrSize = 2,
    ErrRange = 3,
    ErrState = 4,
    ErrUnsupported = 5,
}

#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PhBoundaryAction {
    Reflect = 0,
    Absorb = 1,
}

#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct PhTraceConfig {
    pub max_range_m: PhReal,
    pub max_time_s: PhReal,
    pub surface_depth_m: PhReal,
    pub bottom_depth_m: PhReal,
    pub surface: PhBoundaryAction,
    pub bottom: PhBoundaryAction,
    pub max_bounces: u32,
}

#[repr(C)]
#[derive(Debug, Clone, Copy, Default)]
pub struct PhRayPoint {
    pub range_m: PhReal,
    pub depth_m: PhReal,
    pub angle_rad: PhReal,
    pub time_s: PhReal,
    pub speed_mps: PhReal,
}

#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PhTraceStatus {
    MaxRange = 0,
    MaxTime = 1,
    BufferFull = 2,
    Absorbed = 3,
    MaxBounces = 4,
    Degenerate = 5,
}

#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct PhTraceResult {
    pub point_count: usize,
    pub surface_bounces: u32,
    pub bottom_bounces: u32,
    pub turning_points: u32,
    pub steps: u32,
    pub final_range_m: PhReal,
    pub final_time_s: PhReal,
    pub arc_length_m: PhReal,
    pub status: PhTraceStatus,
}

#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PhPulseType {
    Unknown = 0,
    Cw = 1,
    LfmUp = 2,
    LfmDown = 3,
    Hfm = 4,
}

#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PhTaper {
    Rectangular = 0,
    Hann = 1,
    Tukey25 = 2,
}

#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct PhPulseSpec {
    pub pulse_type: PhPulseType,
    pub f_start_hz: PhReal,
    pub f_end_hz: PhReal,
    pub duration_s: PhReal,
    pub amplitude: PhReal,
    pub taper: PhTaper,
}

#[repr(C)]
#[derive(Debug, Clone, Copy, Default)]
pub struct PhMeasurement {
    pub range_m: PhReal,
    pub bearing_rad: PhReal,
    pub time_s: PhReal,
    pub range_rate_mps: PhReal,
    pub has_range_rate: c_int,
    pub label: u16,
    pub amplitude: PhReal,
}

#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct PhTrackerConfig {
    pub process_accel_mps2: PhReal,
    pub range_sigma_m: PhReal,
    pub bearing_sigma_rad: PhReal,
    pub range_rate_sigma_mps: PhReal,
    pub gate_chi2_2dof: PhReal,
    pub gate_chi2_3dof: PhReal,
    pub confirm_hits: u32,
    pub delete_misses: u32,
    pub init_velocity_sigma_mps: PhReal,
}

#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PhTrackStatus {
    Free = 0,
    Tentative = 1,
    Confirmed = 2,
    Coasting = 3,
}

#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct PhTrackState {
    pub x: PhReal,
    pub y: PhReal,
    pub vx: PhReal,
    pub vy: PhReal,
    pub range_m: PhReal,
    pub bearing_rad: PhReal,
    pub range_rate_mps: PhReal,
    pub id: u32,
    pub hits: u32,
    pub misses: u32,
    pub status: PhTrackStatus,
}

#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PhRsResult {
    Clean = 0,
    Corrected = 1,
    Uncorrectable = 2,
    BadInput = 3,
}

pub const PH_RS_N: usize = 15;
pub const PH_RS_K: usize = 11;

#[repr(C)]
pub struct ph_profile {
    _private: [u8; 0],
}

#[repr(C)]
pub struct ph_tracker {
    _private: [u8; 0],
}

extern "C" {
    pub fn ph_version() -> c_int;
    pub fn ph_version_string() -> *const c_char;
    pub fn ph_real_is_double() -> c_int;
    pub fn ph_status_string(s: PhStatus) -> *const c_char;

    pub fn ph_sound_speed_medwin(t: PhReal, s: PhReal, z: PhReal) -> PhReal;
    pub fn ph_sound_speed_mackenzie(t: PhReal, s: PhReal, z: PhReal) -> PhReal;
    pub fn ph_sound_speed_chen_millero_1977(t: PhReal, s: PhReal, p: PhReal) -> PhReal;
    pub fn ph_sound_speed_chen_millero_its90(t: PhReal, s: PhReal, p: PhReal) -> PhReal;
    pub fn ph_sound_speed_del_grosso(t: PhReal, s: PhReal, p: PhReal) -> PhReal;
    pub fn ph_sound_speed_unesco(t: PhReal, s: PhReal, z: PhReal, lat: PhReal) -> PhReal;
    pub fn ph_depth_to_pressure_bar(z: PhReal, lat: PhReal) -> PhReal;
    pub fn ph_bar_to_kgcm2(bar: PhReal) -> PhReal;
    pub fn ph_t90_to_t68(t90: PhReal) -> PhReal;
    pub fn ph_sound_speed_munk(
        z: PhReal,
        axis_z: PhReal,
        axis_c: PhReal,
        eps: PhReal,
        scale: PhReal,
    ) -> PhReal;

    pub fn ph_profile_size() -> usize;
    pub fn ph_profile_align() -> usize;
    pub fn ph_profile_capacity() -> usize;
    pub fn ph_profile_init(storage: *mut core::ffi::c_void, bytes: usize) -> *mut ph_profile;
    pub fn ph_profile_push(p: *mut ph_profile, depth_m: PhReal, speed_mps: PhReal) -> PhStatus;
    pub fn ph_profile_clear(p: *mut ph_profile);
    pub fn ph_profile_count(p: *const ph_profile) -> usize;
    pub fn ph_profile_sample(
        p: *const ph_profile,
        index: usize,
        out_depth_m: *mut PhReal,
        out_speed_mps: *mut PhReal,
    ) -> PhStatus;
    pub fn ph_profile_speed_at(p: *const ph_profile, depth_m: PhReal) -> PhReal;

    pub fn ph_trace_config_defaults(cfg: *mut PhTraceConfig);
    pub fn ph_trace_ray(
        p: *const ph_profile,
        source_depth_m: PhReal,
        launch_angle_rad: PhReal,
        cfg: *const PhTraceConfig,
        out: *mut PhRayPoint,
        out_capacity: usize,
        result: *mut PhTraceResult,
    ) -> PhStatus;

    pub fn ph_pulse_length(spec: *const PhPulseSpec, fs: PhReal) -> usize;
    pub fn ph_render_real(
        spec: *const PhPulseSpec,
        fs: PhReal,
        out: *mut PhReal,
        out_capacity: usize,
        out_written: *mut usize,
    ) -> PhStatus;
    pub fn ph_render_real_doppler(
        spec: *const PhPulseSpec,
        fs: PhReal,
        doppler: PhReal,
        out: *mut PhReal,
        out_capacity: usize,
        out_written: *mut usize,
    ) -> PhStatus;

    pub fn ph_tracker_config_defaults(cfg: *mut PhTrackerConfig);
    pub fn ph_tracker_size(max_tracks: usize) -> usize;
    pub fn ph_tracker_align() -> usize;
    pub fn ph_tracker_max_tracks() -> usize;
    pub fn ph_tracker_init(
        storage: *mut core::ffi::c_void,
        bytes: usize,
        max_tracks: usize,
    ) -> *mut ph_tracker;
    pub fn ph_tracker_step(
        t: *mut ph_tracker,
        measurements: *const PhMeasurement,
        n: usize,
        cfg: *const PhTrackerConfig,
        time_s: PhReal,
        out_live: *mut usize,
    ) -> PhStatus;
    pub fn ph_tracker_get(t: *const ph_tracker, index: usize, out: *mut PhTrackState) -> PhStatus;
    pub fn ph_tracker_established(t: *const ph_tracker) -> usize;
    pub fn ph_chi2_gate(probability: PhReal, dof: usize) -> PhReal;

    pub fn ph_msequence_length(degree: u32) -> usize;
    pub fn ph_generate_msequence(
        degree: u32,
        seed: u32,
        out: *mut PhReal,
        out_capacity: usize,
        out_written: *mut usize,
    ) -> PhStatus;
    pub fn ph_processing_gain_db(chips_per_bit: usize) -> PhReal;
    pub fn ph_chip_slip(chips_per_bit: usize, closing_speed_mps: PhReal, c_mps: PhReal) -> PhReal;
    pub fn ph_crc32(data: *const u8, length: usize) -> u32;
    pub fn ph_rs_encode(data11: *const u8, out15: *mut u8) -> PhStatus;
    pub fn ph_rs_decode(codeword15: *mut u8, out_corrected: *mut usize) -> PhRsResult;
}
