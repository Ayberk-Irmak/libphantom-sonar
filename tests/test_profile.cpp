// SPDX-License-Identifier: Apache-2.0
#include "framework.hpp"

#include "phantom/profile.hpp"
#include "phantom/sound_speed.hpp"

using namespace phantom;

PT_TEST(profile_rejects_bad_input) {
    SoundSpeedProfile<8> p;
    PT_CHECK(p.push(0, 1500));
    PT_CHECK(!p.push(0, 1510));      // depth must strictly increase
    PT_CHECK(!p.push(-10, 1510));    // ... in the right direction
    PT_CHECK(!p.push(100, -5));      // non-physical sound speed
    PT_CHECK(p.push(100, 1490));
    PT_CHECK(p.valid());
    PT_CHECK(p.size() == 2);
}

PT_TEST(profile_capacity_is_hard) {
    SoundSpeedProfile<4> p;
    for (int i = 0; i < 4; ++i) {
        PT_CHECK(p.push(static_cast<Real>(i * 10), static_cast<Real>(1500 + i)));
    }
    PT_CHECK(!p.push(40, 1504));  // full: refuses rather than allocating
    PT_CHECK(p.size() == 4);
}

PT_TEST(profile_gradients_and_interpolation) {
    SoundSpeedProfile<4> p;
    p.push(0, 1500);
    p.push(100, 1520);   // g = +0.20
    p.push(300, 1500);   // g = -0.10
    const ProfileView v = p.view();

    PT_CHECK(v.valid());
    PT_CHECK(v.layer_count() == 2);
    PT_CHECK_NEAR(v.gradient[0], 0.20, pt::tol(1e-12, 1e-7));
    PT_CHECK_NEAR(v.gradient[1], -0.10, pt::tol(1e-12, 1e-7));

    // Exact at the nodes.
    PT_CHECK_NEAR(speed_at(v, 0),   1500, 1e-12);
    PT_CHECK_NEAR(speed_at(v, 100), 1520, 1e-12);
    PT_CHECK_NEAR(speed_at(v, 300), 1500, 1e-12);
    // Linear in between.
    PT_CHECK_NEAR(speed_at(v, 50),  1510, 1e-12);
    PT_CHECK_NEAR(speed_at(v, 200), 1510, 1e-12);
    // Clamped outside.
    PT_CHECK_NEAR(speed_at(v, -50), 1500, 1e-12);
    PT_CHECK_NEAR(speed_at(v, 900), 1500, 1e-12);
}

PT_TEST(profile_layer_lookup_direction) {
    SoundSpeedProfile<4> p;
    p.push(0, 1500);
    p.push(100, 1520);
    p.push(300, 1500);
    const ProfileView v = p.view();

    // Interior points do not care about direction.
    PT_CHECK(find_layer(v, 50, true) == 0);
    PT_CHECK(find_layer(v, 50, false) == 0);
    PT_CHECK(find_layer(v, 200, true) == 1);
    PT_CHECK(find_layer(v, 200, false) == 1);
    // Sitting exactly on the shared node, direction decides. Getting this wrong
    // is the classic way a ray tracer stalls at a layer interface.
    PT_CHECK(find_layer(v, 100, true) == 1);
    PT_CHECK(find_layer(v, 100, false) == 0);
    // Ends clamp into a real layer.
    PT_CHECK(find_layer(v, 0, false) == 0);
    PT_CHECK(find_layer(v, 300, true) == 1);
}

PT_TEST(profile_fill_from_analytic_function) {
    SoundSpeedProfile<512> p;
    PT_CHECK(fill_profile(p, 0, 5000, 501, [](Real z) { return sound_speed::munk(z); }));
    PT_CHECK(p.size() == 501);
    const ProfileView v = p.view();
    PT_CHECK_NEAR(v.min_depth(), 0, 1e-12);
    PT_CHECK_NEAR(v.max_depth(), 5000, 1e-9);
    PT_CHECK_NEAR(speed_at(v, 1300), 1500.0, 1e-9);
    // Too many samples for the capacity is a refusal, not a truncation.
    SoundSpeedProfile<16> small;
    PT_CHECK(!fill_profile(small, 0, 100, 64, [](Real) { return static_cast<Real>(1500); }));
}
