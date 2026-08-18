/* SPDX-License-Identifier: Apache-2.0
 *
 * The C ABI, exercised from actual C.
 *
 * Compiled with a C compiler in C11 mode, not C++. That is the point: a header
 * that only ever sees a C++ compiler will accumulate C++-isms -- a default
 * argument, a bool, an enum class -- and nobody notices until a real C caller
 * arrives. This file is the thing that notices.
 */
#include "phantom/phantom.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int g_failures = 0;
static int g_checks = 0;

#define CHECK(cond)                                                            \
    do {                                                                       \
        ++g_checks;                                                            \
        if (!(cond)) {                                                         \
            printf("    FAIL %s:%d  %s\n", __FILE__, __LINE__, #cond);         \
            ++g_failures;                                                      \
        }                                                                      \
    } while (0)

#define CHECK_NEAR(got, want, tol)                                             \
    do {                                                                       \
        ++g_checks;                                                            \
        if (!(fabs((double)(got) - (double)(want)) <= (tol))) {                \
            printf("    FAIL %s:%d  %s ~= %s (got %.6f, want %.6f)\n",         \
                   __FILE__, __LINE__, #got, #want, (double)(got),             \
                   (double)(want));                                            \
            ++g_failures;                                                      \
        }                                                                      \
    } while (0)

static void test_version_and_precision(void) {
    printf("  version %s (%d), Real is %s\n", ph_version_string(), ph_version(),
           ph_real_is_double() ? "double" : "float");
    CHECK(ph_version() >= 1400);
    /* The check every C caller should make at startup: a library built with a
     * different Real than the caller assumed reinterprets every number crossing
     * the boundary, which produces plausible nonsense rather than a crash. */
    CHECK(ph_real_is_double() == PH_REAL_IS_DOUBLE);
    CHECK(strcmp(ph_status_string(PH_OK), "ok") == 0);
}

static void test_sound_speed(void) {
    /* The published UNESCO check value, through the C boundary. */
    ph_real c = ph_sound_speed_chen_millero_1977((ph_real)40, (ph_real)40, (ph_real)1000);
    printf("  UNESCO 44 check value through the C ABI: %.4f (published 1731.995)\n", (double)c);
    CHECK_NEAR(c, 1731.995, PH_REAL_IS_DOUBLE ? 1e-3 : 0.05);

    /* The three equations must agree where they overlap. */
    ph_real med = ph_sound_speed_medwin((ph_real)10, (ph_real)35, (ph_real)500);
    ph_real mac = ph_sound_speed_mackenzie((ph_real)10, (ph_real)35, (ph_real)500);
    CHECK_NEAR(med, mac, 1.0);

    ph_real bar = ph_depth_to_pressure_bar((ph_real)1000, (ph_real)45);
    CHECK(bar > 100.0 && bar < 102.0);
}

static void test_profile_and_ray(void) {
    /* Caller-owned storage: this is the pattern the whole ABI is built on. */
    size_t bytes = ph_profile_size();
    void* storage = malloc(bytes);
    CHECK(storage != NULL);
    ph_profile* p = ph_profile_init(storage, bytes);
    CHECK(p != NULL);
    CHECK(ph_profile_init(storage, bytes - 1) == NULL);   /* too small, refused */
    CHECK(ph_profile_init(NULL, bytes) == NULL);

    for (int i = 0; i <= 200; ++i) {
        ph_real z = (ph_real)(i * 25);
        ph_real c = ph_sound_speed_munk(z, (ph_real)1300, (ph_real)1500,
                                        (ph_real)7.37e-3, (ph_real)1300);
        CHECK(ph_profile_push(p, z, c) == PH_OK);
    }
    CHECK(ph_profile_count(p) == 201);
    /* Depths must strictly increase; going backwards is refused, not accepted. */
    CHECK(ph_profile_push(p, (ph_real)0, (ph_real)1500) == PH_ERR_RANGE);

    ph_real z0 = 0, c0 = 0;
    CHECK(ph_profile_sample(p, 52, &z0, &c0) == PH_OK);
    CHECK_NEAR(z0, 1300.0, 1e-6);
    CHECK_NEAR(c0, 1500.0, 1e-6);
    CHECK(ph_profile_sample(p, 999, &z0, &c0) == PH_ERR_RANGE);

    ph_trace_config cfg;
    ph_trace_config_defaults(&cfg);
    cfg.max_range_m = (ph_real)40000;
    cfg.max_time_s = (ph_real)40;

    static ph_ray_point points[4096];
    ph_trace_result r;
    CHECK(ph_trace_ray(p, (ph_real)1300, (ph_real)0.12, &cfg, points, 4096, &r) == PH_OK);
    printf("  ray from the axis: %zu points, %.0f m, %u turning points\n",
           r.point_count, (double)r.final_range_m, r.turning_points);
    CHECK(r.point_count > 10);
    CHECK(r.turning_points > 0);

    /* Snell's invariant across the traced path. */
    double xi0 = cos((double)points[0].angle_rad) / (double)points[0].speed_mps;
    double worst = 0;
    for (size_t i = 1; i < r.point_count; ++i) {
        double xi = cos((double)points[i].angle_rad) / (double)points[i].speed_mps;
        double d = fabs((xi - xi0) / xi0);
        if (d > worst) worst = d;
    }
    printf("  worst Snell drift through the C ABI: %.2e\n", worst);
    CHECK(worst < (PH_REAL_IS_DOUBLE ? 1e-12 : 1e-5));

    free(storage);
}

static void test_waveform(void) {
    ph_pulse_spec spec;
    spec.type = PH_PULSE_HFM;
    spec.f_start_hz = (ph_real)18000;
    spec.f_end_hz = (ph_real)30000;
    spec.duration_s = (ph_real)0.01;
    spec.amplitude = (ph_real)1;
    spec.taper = PH_TAPER_TUKEY25;

    size_t n = ph_pulse_length(&spec, (ph_real)96000);
    CHECK(n == 960);

    static ph_real samples[4096];
    size_t written = 0;
    CHECK(ph_render_real(&spec, (ph_real)96000, samples, 4096, &written) == PH_OK);
    CHECK(written == 960);
    /* A tapered pulse starts near zero and has energy in the middle. */
    CHECK(fabs((double)samples[0]) < 0.2);
    double energy = 0;
    for (size_t i = 0; i < written; ++i) energy += (double)samples[i] * (double)samples[i];
    CHECK(energy > 100.0);

    /* Too small a buffer must be refused, not overrun. */
    CHECK(ph_render_real(&spec, (ph_real)96000, samples, 10, &written) == PH_ERR_SIZE);
}

static void test_tracker(void) {
    size_t bytes = ph_tracker_size(8);
    CHECK(bytes > 0);
    CHECK(ph_tracker_size(0) == 0);
    CHECK(ph_tracker_size(ph_tracker_max_tracks() + 1) == 0);

    void* storage = malloc(bytes);
    ph_tracker* t = ph_tracker_init(storage, bytes, 8);
    CHECK(t != NULL);

    ph_tracker_config cfg;
    ph_tracker_config_defaults(&cfg);
    cfg.range_sigma_m = (ph_real)5;
    /* The target below closes at 54 m/s and the default initial velocity sigma
     * is 10, which would pin the estimate near zero for many scans. A caller
     * must size this to the targets it expects -- the first version of this test
     * did not, and the tracker converged to half the true closing rate while
     * reporting a perfectly healthy track. */
    cfg.init_velocity_sigma_mps = (ph_real)60;

    /* A target closing straight in. */
    double x = 1500, y = 7000, vx = -20, vy = -50;
    for (int k = 1; k <= 20; ++k) {
        x += vx; y += vy;
        ph_measurement m;
        memset(&m, 0, sizeof m);
        m.range_m = (ph_real)sqrt(x * x + y * y);
        m.bearing_rad = (ph_real)atan2(x, y);
        m.time_s = (ph_real)k;
        m.amplitude = (ph_real)1;
        size_t live = 0;
        CHECK(ph_tracker_step(t, &m, 1, &cfg, (ph_real)k, &live) == PH_OK);
    }
    size_t established = ph_tracker_established(t);
    printf("  20 clean detections -> %zu established track(s)\n", established);
    CHECK(established == 1);

    for (size_t i = 0; i < 8; ++i) {
        ph_track_state s;
        if (ph_tracker_get(t, i, &s) != PH_OK) continue;
        if (s.status != PH_TRACK_CONFIRMED && s.status != PH_TRACK_COASTING) continue;
        printf("  track %u: range %.0f m, closing %.1f m/s (truth %.1f)\n",
               s.id, (double)s.range_m, (double)s.range_rate_mps,
               -(x * vx + y * vy) / sqrt(x * x + y * y));
        CHECK_NEAR(s.range_rate_mps, -(x * vx + y * vy) / sqrt(x * x + y * y), 6.0);
    }
    CHECK(ph_tracker_get(t, 99, NULL) == PH_ERR_NULL);

    CHECK_NEAR(ph_chi2_gate((ph_real)0.95, 2), 5.991, 1e-2);
    CHECK_NEAR(ph_chi2_gate((ph_real)0.95, 3), 7.815, 1e-2);
    free(storage);
}

static void test_comm(void) {
    CHECK(ph_msequence_length(9) == 511);
    CHECK(ph_msequence_length(3) == 0);

    static ph_real chips[1024];
    size_t written = 0;
    CHECK(ph_generate_msequence(9, 0x1FFu, chips, 1024, &written) == PH_OK);
    CHECK(written == 511);
    /* Balance: exactly one more +1 than -1 over a full period. */
    long sum = 0;
    for (size_t i = 0; i < written; ++i) sum += (chips[i] > 0) ? 1 : -1;
    printf("  511-chip m-sequence balance through the C ABI: %ld\n", sum);
    CHECK(sum == 1 || sum == -1);
    CHECK(ph_generate_msequence(9, 0u, chips, 1024, &written) == PH_ERR_SIZE);
    CHECK(ph_generate_msequence(3, 1u, chips, 1024, &written) == PH_ERR_RANGE);

    /* The published CRC-32 check value. */
    const unsigned char s[9] = {'1','2','3','4','5','6','7','8','9'};
    unsigned int crc = ph_crc32(s, 9);
    printf("  CRC-32(\"123456789\") through the C ABI = 0x%08X\n", crc);
    CHECK(crc == 0xCBF43926u);

    /* Reed-Solomon round trip with two symbol errors. */
    unsigned char data[PH_RS_K];
    unsigned char code[PH_RS_N];
    for (int i = 0; i < PH_RS_K; ++i) data[i] = (unsigned char)(i % 16);
    CHECK(ph_rs_encode(data, code) == PH_OK);
    unsigned char rx[PH_RS_N];
    memcpy(rx, code, PH_RS_N);
    rx[2] ^= 0x0B;
    rx[9] ^= 0x07;
    size_t fixed = 0;
    ph_rs_result r = ph_rs_decode(rx, &fixed);
    printf("  RS(15,11) with 2 symbol errors -> %s, %zu repaired\n",
           (r == PH_RS_CORRECTED) ? "corrected" : "failed", fixed);
    CHECK(r == PH_RS_CORRECTED);
    CHECK(fixed == 2);
    CHECK(memcmp(rx, code, PH_RS_N) == 0);

    printf("  chip slip, 511 chips at 3 m/s: %.3f chips\n",
           (double)ph_chip_slip(511, (ph_real)3, (ph_real)1500));
    CHECK(ph_chip_slip(511, (ph_real)3, (ph_real)1500) > 1.0);
}

int main(void) {
    printf("libphantom-sonar C ABI\n");
    printf("--------------------------------------------------------------\n");
    test_version_and_precision();
    test_sound_speed();
    test_profile_and_ray();
    test_waveform();
    test_tracker();
    test_comm();
    printf("--------------------------------------------------------------\n");
    printf("%d checks, %d failures\n", g_checks, g_failures);
    return g_failures == 0 ? 0 : 1;
}
