// SPDX-License-Identifier: Apache-2.0
//! Safe Rust bindings for `libphantom-sonar`.
//!
//! # What "safe" buys here, and what it does not
//!
//! The C library allocates nothing and every stateful object lives in storage
//! the caller supplies. Rust makes that ownership expressible rather than
//! documented: a [`Profile`] owns its own boxed storage, its lifetime is the
//! compiler's problem, and there is no way to hand the library a buffer that
//! outlives it or is too small.
//!
//! What these bindings do NOT do is re-verify the physics. Every number here
//! comes from the same C++ code the test suite measures, and the tests in this
//! crate check the *binding* -- that a value survives the round trip, that a
//! size mismatch is refused, that a struct's layout matches. A ray path is not
//! more correct for having been fetched through Rust.
//!
//! # Precision
//!
//! Build with the `real_float` feature if and only if the C library was built
//! with `PHANTOM_REAL_FLOAT`. A mismatch is not a link error; it reinterprets
//! every float crossing the boundary. [`check_precision`] compares the two at
//! runtime and every constructor calls it, so the failure is loud instead of
//! silent.
use phantom_sonar_sys as sys;
use std::ffi::CStr;

pub use sys::PhReal as Real;

/// Errors the C ABI can report.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Error {
    /// A required pointer was null.
    Null,
    /// A buffer was too small, or a size was zero.
    Size,
    /// An argument was outside its valid range.
    Range,
    /// The object was not in a state allowing the call.
    State,
    /// The build does not provide this.
    Unsupported,
    /// The library was built with a different `Real` than this crate assumes.
    /// Every float crossing the boundary would be reinterpreted.
    PrecisionMismatch { library_is_double: bool, crate_is_double: bool },
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Error::PrecisionMismatch { library_is_double, crate_is_double } => write!(
                f,
                "precision mismatch: library Real is {}, crate assumes {}. \
                 Rebuild one to match; the `real_float` cargo feature selects f32.",
                if *library_is_double { "double" } else { "float" },
                if *crate_is_double { "f64" } else { "f32" },
            ),
            other => write!(f, "{other:?}"),
        }
    }
}

impl std::error::Error for Error {}

pub type Result<T> = std::result::Result<T, Error>;

fn check(status: sys::PhStatus) -> Result<()> {
    match status {
        sys::PhStatus::Ok => Ok(()),
        sys::PhStatus::ErrNull => Err(Error::Null),
        sys::PhStatus::ErrSize => Err(Error::Size),
        sys::PhStatus::ErrRange => Err(Error::Range),
        sys::PhStatus::ErrState => Err(Error::State),
        sys::PhStatus::ErrUnsupported => Err(Error::Unsupported),
    }
}

/// `true` when this crate's `Real` matches the linked library's.
pub fn check_precision() -> Result<()> {
    let library_is_double = unsafe { sys::ph_real_is_double() } != 0;
    let crate_is_double = std::mem::size_of::<Real>() == 8;
    if library_is_double == crate_is_double {
        Ok(())
    } else {
        Err(Error::PrecisionMismatch { library_is_double, crate_is_double })
    }
}

/// Library version as `major * 10000 + minor * 100 + patch`.
pub fn version() -> i32 {
    unsafe { sys::ph_version() }
}

/// Library version as a string.
pub fn version_string() -> &'static str {
    unsafe { CStr::from_ptr(sys::ph_version_string()) }.to_str().unwrap_or("?")
}

/// Sound speed in seawater. See `docs/math_spec.md` §1 and §17 for which
/// equation to use; they are not interchangeable and they disagree by more than
/// their individual stated uncertainties.
pub mod sound_speed {
    use super::Real;
    use phantom_sonar_sys as sys;

    /// Medwin (1975). Valid to about 1000 m, which is at or above the deep
    /// sound channel axis at most latitudes -- so this is the wrong tool for
    /// long-range work.
    pub fn medwin(temperature_c: Real, salinity_psu: Real, depth_m: Real) -> Real {
        unsafe { sys::ph_sound_speed_medwin(temperature_c, salinity_psu, depth_m) }
    }

    /// Mackenzie (1981), nine-term. Valid to 8000 m.
    pub fn mackenzie(temperature_c: Real, salinity_psu: Real, depth_m: Real) -> Real {
        unsafe { sys::ph_sound_speed_mackenzie(temperature_c, salinity_psu, depth_m) }
    }

    /// Original Chen & Millero (1977), IPTS-68, pressure in bar. This is the
    /// version the published UNESCO check value and table belong to.
    pub fn chen_millero_1977(temperature_c: Real, salinity_psu: Real, pressure_bar: Real) -> Real {
        unsafe { sys::ph_sound_speed_chen_millero_1977(temperature_c, salinity_psu, pressure_bar) }
    }

    /// Wong & Zhu (1995), ITS-90. What modern data wants.
    pub fn chen_millero_its90(temperature_c: Real, salinity_psu: Real, pressure_bar: Real) -> Real {
        unsafe { sys::ph_sound_speed_chen_millero_its90(temperature_c, salinity_psu, pressure_bar) }
    }

    /// Del Grosso (1974), ITS-90 form, pressure in kg/cm².
    pub fn del_grosso(temperature_c: Real, salinity_psu: Real, pressure_kgcm2: Real) -> Real {
        unsafe { sys::ph_sound_speed_del_grosso(temperature_c, salinity_psu, pressure_kgcm2) }
    }

    /// UNESCO from depth rather than pressure.
    pub fn unesco(temperature_c: Real, salinity_psu: Real, depth_m: Real, latitude_deg: Real) -> Real {
        unsafe { sys::ph_sound_speed_unesco(temperature_c, salinity_psu, depth_m, latitude_deg) }
    }

    /// Leroy & Parthiot depth to gauge pressure, in bar.
    pub fn depth_to_pressure_bar(depth_m: Real, latitude_deg: Real) -> Real {
        unsafe { sys::ph_depth_to_pressure_bar(depth_m, latitude_deg) }
    }

    /// Munk's canonical deep sound channel.
    pub fn munk(depth_m: Real, axis_depth_m: Real, axis_speed_mps: Real, epsilon: Real, scale_m: Real) -> Real {
        unsafe { sys::ph_sound_speed_munk(depth_m, axis_depth_m, axis_speed_mps, epsilon, scale_m) }
    }
}

/// A sound speed profile, owning its storage.
pub struct Profile {
    // Never read, and that is the point: this field exists to OWN the bytes
    // `handle` points into. Dropping it frees them, and nothing else may. This
    // is the Rust expression of the C ABI's caller-owned-storage contract --
    // the C side has no destructor because it allocated nothing.
    #[allow(dead_code)]
    storage: Box<[u8]>,
    handle: *mut sys::ph_profile,
}

// The handle points into `storage`, which this struct owns exclusively, and the
// C side holds no global state. Sending one to another thread moves the storage
// with it.
unsafe impl Send for Profile {}

impl Profile {
    /// Allocates storage and places a profile in it.
    pub fn new() -> Result<Self> {
        check_precision()?;
        let bytes = unsafe { sys::ph_profile_size() };
        let align = unsafe { sys::ph_profile_align() };
        // Over-allocate and align by hand: Box<[u8]> is only byte-aligned, and
        // the C side refuses a misaligned pointer rather than accepting one and
        // producing an unaligned load.
        let mut storage = vec![0u8; bytes + align].into_boxed_slice();
        let base = storage.as_mut_ptr() as usize;
        let offset = (align - (base % align)) % align;
        let handle = unsafe {
            sys::ph_profile_init(storage.as_mut_ptr().add(offset) as *mut _, bytes)
        };
        if handle.is_null() {
            return Err(Error::Size);
        }
        Ok(Profile { storage, handle })
    }

    /// Maximum samples a profile can hold. Fixed at library build time.
    pub fn capacity() -> usize {
        unsafe { sys::ph_profile_capacity() }
    }

    /// Appends a sample. Depths must strictly increase.
    pub fn push(&mut self, depth_m: Real, speed_mps: Real) -> Result<()> {
        check(unsafe { sys::ph_profile_push(self.handle, depth_m, speed_mps) })
    }

    pub fn clear(&mut self) {
        unsafe { sys::ph_profile_clear(self.handle) }
    }

    pub fn len(&self) -> usize {
        unsafe { sys::ph_profile_count(self.handle) }
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// `(depth, speed)` at an index.
    pub fn sample(&self, index: usize) -> Result<(Real, Real)> {
        let mut z: Real = 0.0;
        let mut c: Real = 0.0;
        check(unsafe { sys::ph_profile_sample(self.handle, index, &mut z, &mut c) })?;
        Ok((z, c))
    }

    /// Linear interpolation, clamped outside the profile.
    ///
    /// A profile needs at least two samples to interpolate between. The C entry
    /// point returns 0 m/s for a shorter one, which is a silent failure a caller
    /// can easily propagate; the safe layer turns it into an error instead.
    /// Catching this is the sort of thing the wrapper is for.
    pub fn speed_at(&self, depth_m: Real) -> Result<Real> {
        if self.len() < 2 {
            return Err(Error::State);
        }
        Ok(unsafe { sys::ph_profile_speed_at(self.handle, depth_m) })
    }

    /// Traces one ray, returning its path and a summary.
    pub fn trace_ray(
        &self,
        source_depth_m: Real,
        launch_angle_rad: Real,
        config: &TraceConfig,
        max_points: usize,
    ) -> Result<(Vec<RayPoint>, TraceSummary)> {
        let mut points: Vec<sys::PhRayPoint> = vec![Default::default(); max_points];
        let mut result = std::mem::MaybeUninit::<sys::PhTraceResult>::uninit();
        check(unsafe {
            sys::ph_trace_ray(
                self.handle,
                source_depth_m,
                launch_angle_rad,
                &config.0,
                points.as_mut_ptr(),
                max_points,
                result.as_mut_ptr(),
            )
        })?;
        let result = unsafe { result.assume_init() };
        points.truncate(result.point_count);
        let path = points
            .into_iter()
            .map(|p| RayPoint {
                range_m: p.range_m,
                depth_m: p.depth_m,
                angle_rad: p.angle_rad,
                time_s: p.time_s,
                speed_mps: p.speed_mps,
            })
            .collect();
        Ok((
            path,
            TraceSummary {
                surface_bounces: result.surface_bounces,
                bottom_bounces: result.bottom_bounces,
                turning_points: result.turning_points,
                final_range_m: result.final_range_m,
                final_time_s: result.final_time_s,
                arc_length_m: result.arc_length_m,
            },
        ))
    }
}

/// One point on a traced ray.
#[derive(Debug, Clone, Copy)]
pub struct RayPoint {
    pub range_m: Real,
    pub depth_m: Real,
    /// Grazing angle, positive downgoing.
    pub angle_rad: Real,
    pub time_s: Real,
    pub speed_mps: Real,
}

/// What a trace did, apart from where it went.
#[derive(Debug, Clone, Copy)]
pub struct TraceSummary {
    pub surface_bounces: u32,
    pub bottom_bounces: u32,
    pub turning_points: u32,
    pub final_range_m: Real,
    pub final_time_s: Real,
    pub arc_length_m: Real,
}

/// Ray tracing limits and boundary behaviour.
pub struct TraceConfig(sys::PhTraceConfig);

impl Default for TraceConfig {
    fn default() -> Self {
        let mut cfg = std::mem::MaybeUninit::<sys::PhTraceConfig>::uninit();
        unsafe { sys::ph_trace_config_defaults(cfg.as_mut_ptr()) };
        TraceConfig(unsafe { cfg.assume_init() })
    }
}

impl TraceConfig {
    pub fn max_range_m(mut self, v: Real) -> Self {
        self.0.max_range_m = v;
        self
    }
    pub fn max_time_s(mut self, v: Real) -> Self {
        self.0.max_time_s = v;
        self
    }
    pub fn absorb_at_boundaries(mut self) -> Self {
        self.0.surface = sys::PhBoundaryAction::Absorb;
        self.0.bottom = sys::PhBoundaryAction::Absorb;
        self
    }
}

/// Error detection and correction for an acoustic link.
pub mod coding {
    use phantom_sonar_sys as sys;

    /// CRC-32 (IEEE 802.3, reflected).
    pub fn crc32(data: &[u8]) -> u32 {
        unsafe { sys::ph_crc32(data.as_ptr(), data.len()) }
    }

    /// Reed-Solomon (15, 11) over GF(16). Corrects up to 2 symbol errors.
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub enum RsResult {
        /// No errors detected.
        Clean,
        /// Errors found and repaired.
        Corrected(usize),
        /// Errors detected, too many to repair. This is a real answer: a
        /// decoder that always returns something produces plausible wrong data.
        Uncorrectable,
    }

    /// Encodes 11 data symbols (each < 16) into a 15-symbol codeword.
    pub fn rs_encode(data: &[u8; sys::PH_RS_K]) -> Option<[u8; sys::PH_RS_N]> {
        let mut out = [0u8; sys::PH_RS_N];
        let ok = unsafe { sys::ph_rs_encode(data.as_ptr(), out.as_mut_ptr()) };
        (ok == sys::PhStatus::Ok).then_some(out)
    }

    /// Decodes in place.
    pub fn rs_decode(codeword: &mut [u8; sys::PH_RS_N]) -> RsResult {
        let mut fixed = 0usize;
        match unsafe { sys::ph_rs_decode(codeword.as_mut_ptr(), &mut fixed) } {
            sys::PhRsResult::Clean => RsResult::Clean,
            sys::PhRsResult::Corrected => RsResult::Corrected(fixed),
            _ => RsResult::Uncorrectable,
        }
    }
}
