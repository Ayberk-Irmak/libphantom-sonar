// SPDX-License-Identifier: Apache-2.0
//
// These tests check the BINDING, not the physics.
//
// The physics is measured by the C++ suite, against closed forms and published
// tables. Re-asserting a ray path here would prove nothing except that FFI
// copies bytes. What can go wrong at this layer is different and specific:
// a struct field in the wrong order, an enum discriminant off by one, a size or
// alignment mismatch, a precision mismatch that silently reinterprets floats.
// So the checks below either round-trip a value with a known answer, or compare
// a layout against what the C library itself reports.
use phantom_sonar::{coding, sound_speed, Error, Profile, TraceConfig};

#[test]
fn precision_matches_the_linked_library() {
    // The check every binding must make first. A library built with float and a
    // crate assuming f64 is not a link error -- it reinterprets every number
    // crossing the boundary and produces plausible nonsense.
    phantom_sonar::check_precision().expect("build the crate with --features real_float to match");
}

#[test]
fn version_is_reported_through_the_boundary() {
    assert!(phantom_sonar::version() >= 10000);
    assert!(!phantom_sonar::version_string().is_empty());
}

#[test]
fn the_published_unesco_check_value_survives_the_round_trip() {
    // Not a test of Chen-Millero -- that is verified against the primary source
    // in the C++ suite. It is a test that a f64 goes in and comes back with its
    // digits intact, using a value whose correct answer is published.
    let c = sound_speed::chen_millero_1977(40.0, 40.0, 1000.0);
    let tol = if std::mem::size_of::<phantom_sonar::Real>() == 8 { 1e-3 } else { 0.05 };
    assert!(
        (c as f64 - 1731.995).abs() < tol,
        "UNESCO 44 check value came back as {c}, expected 1731.995"
    );
}

#[test]
fn a_profile_owns_its_storage_and_traces_a_ray() {
    let mut p = Profile::new().unwrap();
    assert!(Profile::capacity() >= 201);
    for i in 0..=200 {
        let z = (i * 25) as phantom_sonar::Real;
        let c = sound_speed::munk(z, 1300.0, 1500.0, 7.37e-3, 1300.0);
        p.push(z, c).unwrap();
    }
    assert_eq!(p.len(), 201);

    // Depths must strictly increase; the error is surfaced as a Result rather
    // than a silently ignored sample.
    assert_eq!(p.push(0.0, 1500.0), Err(Error::Range));
    assert_eq!(p.sample(999).unwrap_err(), Error::Range);

    let (z, c) = p.sample(52).unwrap();
    assert!((z as f64 - 1300.0).abs() < 1e-6);
    assert!((c as f64 - 1500.0).abs() < 1e-6);

    let cfg = TraceConfig::default().max_range_m(40_000.0).max_time_s(40.0);
    let (path, summary) = p.trace_ray(1300.0, 0.12, &cfg, 4096).unwrap();
    assert!(path.len() > 10);
    assert!(summary.turning_points > 0);

    // Snell's invariant across the path, computed in Rust from the returned
    // points. This is a layout check as much as a physics one: if the struct
    // fields were transposed, angle and speed would be swapped and this fails.
    let xi0 = (path[0].angle_rad as f64).cos() / path[0].speed_mps as f64;
    let worst = path
        .iter()
        .map(|pt| (((pt.angle_rad as f64).cos() / pt.speed_mps as f64) - xi0).abs() / xi0)
        .fold(0.0f64, f64::max);
    assert!(worst < 1e-5, "Snell drift through the binding: {worst:e}");
}

#[test]
fn a_dropped_profile_does_not_take_a_live_handle_with_it() {
    // Ownership is the whole reason for the safe layer: the C side has no
    // destructor because it allocated nothing, so the storage must die with the
    // Rust object and not before.
    let handle_len = {
        let mut p = Profile::new().unwrap();
        p.push(0.0, 1500.0).unwrap();
        p.push(100.0, 1510.0).unwrap();
        p.len()
    };
    assert_eq!(handle_len, 2);
    // A second profile reuses freed memory; if the first had left anything
    // dangling this is where it would show under a sanitizer.
    let mut q = Profile::new().unwrap();
    q.push(0.0, 1490.0).unwrap();
    assert_eq!(q.len(), 1);
    // One sample is not a profile: there is nothing to interpolate between. The
    // C entry point answers 0 m/s here, which a caller could propagate without
    // noticing; the safe layer refuses instead.
    assert_eq!(q.speed_at(50.0), Err(Error::State));
    q.push(100.0, 1510.0).unwrap();
    assert!((q.speed_at(50.0).unwrap() as f64 - 1500.0).abs() < 1e-9);
}

#[test]
fn crc32_matches_its_published_check_value() {
    assert_eq!(coding::crc32(b"123456789"), 0xCBF4_3926);
    assert_eq!(coding::crc32(b""), 0);
}

#[test]
fn reed_solomon_round_trips_through_the_binding() {
    let data: [u8; 11] = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
    let code = coding::rs_encode(&data).expect("11 symbols under 16 must encode");
    // Systematic: the data is the front of the codeword.
    assert_eq!(&code[..11], &data[..]);

    let mut rx = code;
    rx[2] ^= 0x0B;
    rx[9] ^= 0x07;
    assert_eq!(coding::rs_decode(&mut rx), coding::RsResult::Corrected(2));
    assert_eq!(rx, code);

    // A symbol out of range must be refused rather than silently masked.
    let bad: [u8; 11] = [99, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
    assert!(coding::rs_encode(&bad).is_none());
}

#[test]
fn struct_layouts_agree_with_the_c_library() {
    // The hand-written declarations in the -sys crate can drift from phantom.h.
    // Comparing sizes and alignments against what the library itself reports is
    // what turns drift into a test failure rather than a memory corruption.
    use phantom_sonar_sys as sys;
    unsafe {
        assert_eq!(std::mem::size_of::<sys::ph_profile>(), 0, "opaque types must be zero-sized");
        assert!(sys::ph_profile_size() > 0);
        assert!(sys::ph_profile_align() >= std::mem::align_of::<sys::PhReal>());
        assert!(sys::ph_tracker_size(8) > 0);
        assert_eq!(sys::ph_tracker_size(0), 0, "zero tracks must be refused");
        assert_eq!(
            sys::ph_tracker_size(sys::ph_tracker_max_tracks() + 1),
            0,
            "over-capacity must be refused"
        );
    }
    // Enums must be C-sized, or a struct containing one is the wrong length.
    assert_eq!(std::mem::size_of::<sys::PhTraceStatus>(), 4);
    assert_eq!(std::mem::size_of::<sys::PhStatus>(), 4);
}
