// SPDX-License-Identifier: Apache-2.0
//
// Finds the built C library. No bindgen: the declarations in lib.rs are written
// by hand against phantom.h, for the same reason phantom.h is written by hand
// against the C++ headers -- a generated binding tracks whatever the header
// happens to say today, including its mistakes, and gives no place to record
// what a function actually means.
//
// The cost is that the two can drift. `cargo test` runs a layout test that
// compares every struct's size and alignment against the C library's own view
// of them, which is what turns drift into a failure instead of a corruption.
use std::env;

fn main() {
    let dir = env::var("PHANTOM_SONAR_LIB_DIR")
        .unwrap_or_else(|_| "../../../build".to_string());
    println!("cargo:rustc-link-search=native={dir}");
    println!("cargo:rustc-link-lib=static=phantom-sonar");
    // The library is C++ underneath, so its runtime must come along.
    println!("cargo:rustc-link-lib=dylib=stdc++");
    println!("cargo:rustc-link-lib=dylib=m");
    println!("cargo:rerun-if-env-changed=PHANTOM_SONAR_LIB_DIR");
}
