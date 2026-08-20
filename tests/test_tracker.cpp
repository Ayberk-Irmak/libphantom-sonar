// SPDX-License-Identifier: Apache-2.0
// Detection-to-track association and target state.
//
// The load-bearing test is the NIS consistency check. A Kalman filter that
// tracks beautifully and reports the wrong covariance passes every position
// test ever written and is still broken -- it will gate correct measurements
// out, or accept clutter, and nothing about its output says so. The normalised
// innovation squared has a known distribution when the filter is right, and
// checking it is the standard way to find out.
#include "framework.hpp"

#include "phantom/tracker.hpp"

#include <algorithm>
#include <array>
#include <cmath>

using namespace phantom;

namespace {

constexpr double kPiD = 3.14159265358979323846;

std::array<Track, 32> g_tracks;
std::array<Measurement, 16> g_meas;

TrackerConfig default_cfg() {
    TrackerConfig cfg;
    cfg.process_accel_mps2 = static_cast<Real>(0.2);
    cfg.range_sigma_m = 10;
    cfg.bearing_sigma_rad = deg2rad(static_cast<Real>(1));
    cfg.gate_chi2_2dof = chi2_gate_2dof(static_cast<Real>(0.99));
    cfg.confirm_hits = 3;
    cfg.delete_misses = 3;
    return cfg;
}

void clear_tracks() {
    for (Track& t : g_tracks) t = Track{};
}

// A constant-velocity truth target.
struct Truth {
    double x, y, vx, vy;
    void advance(double dt) { x += vx * dt; y += vy * dt; }
    double range() const { return std::sqrt(x * x + y * y); }
    double bearing() const { return std::atan2(x, y); }
};

Measurement observe(const Truth& t, const TrackerConfig& cfg, pt::Rng& rng, Real time_s) {
    Measurement z;
    z.range_m = static_cast<Real>(t.range()
              + static_cast<double>(cfg.range_sigma_m) * rng.normal());
    z.bearing_rad = static_cast<Real>(t.bearing()
                  + static_cast<double>(cfg.bearing_sigma_rad) * rng.normal());
    z.time_s = time_s;
    return z;
}

}  // namespace

// ---------------------------------------------------------------------------
// Geometry and the gate
// ---------------------------------------------------------------------------

PT_TEST(target_state_converts_between_cartesian_and_polar) {
    TargetState s;
    s.x = 300; s.y = 400; s.vx = -3; s.vy = -4;
    PT_CHECK_REL(s.range_m(), 500.0, pt::tol(1e-12, 1e-5));
    PT_CHECK_NEAR(rad2deg(s.bearing_rad()), std::atan2(300.0, 400.0) * 180.0 / kPiD,
                  pt::tol(1e-9, 1e-4));
    // Moving straight at the origin at 5 m/s: the closing rate is the speed.
    PT_CHECK_REL(s.range_rate_mps(), 5.0, pt::tol(1e-9, 1e-4));

    // Crossing at constant range: no closing at all.
    TargetState c;
    c.x = 0; c.y = 1000; c.vx = 7; c.vy = 0;
    PT_CHECK_NEAR(c.range_rate_mps(), 0.0, pt::tol(1e-9, 1e-3));
    PT_CHECK_NEAR(c.bearing_rad(), 0.0, 1e-12);

    TargetState origin;
    PT_CHECK(origin.range_rate_mps() == 0);
}

PT_TEST(chi_square_gate_inverts_in_closed_form) {
    // The 2-dof chi-square CDF is 1 - exp(-x/2), so the gate needs no table.
    PT_CHECK_NEAR(chi2_gate_2dof(static_cast<Real>(0.95)), 5.991, 0.002);
    PT_CHECK_NEAR(chi2_gate_2dof(static_cast<Real>(0.99)), 9.210, 0.002);
    PT_CHECK_NEAR(chi2_gate_2dof(static_cast<Real>(0.5)), 1.386, 0.002);
    PT_CHECK(chi2_gate_2dof(0) == 0);
    PT_CHECK(chi2_gate_2dof(1) == 0);
    PT_CHECK(chi2_gate_2dof(2) == 0);
}

// ---------------------------------------------------------------------------
// Filter consistency -- the test that matters
// ---------------------------------------------------------------------------

PT_TEST(filter_nis_follows_its_chi_square_distribution) {
    // A filter whose covariance is wrong still tracks; it just lies about how
    // well. NIS is chi-square with 2 dof when the filter is consistent, so its
    // mean must be 2 and 95% of samples must fall under 5.991.
    const TrackerConfig cfg = default_cfg();
    pt::Rng rng(20260807);

    double sum_nis = 0;
    std::size_t n = 0;
    std::size_t under_95 = 0;
    std::size_t under_99 = 0;

    const std::size_t runs = 60;
    for (std::size_t r = 0; r < runs; ++r) {
        Truth truth{ -200.0 + 8.0 * rng.uniform(), 2500.0 + 400.0 * rng.uniform(),
                      3.0 * rng.uniform(), -6.0 + 2.0 * rng.uniform() };
        Track t{};
        track_initiate(t, observe(truth, cfg, rng, 0), cfg, 1);

        const double dt = 1.0;
        for (std::size_t k = 1; k < 40; ++k) {
            truth.advance(dt);
            track_predict(t, static_cast<Real>(dt), cfg);
            const Measurement z = observe(truth, cfg, rng, static_cast<Real>(static_cast<double>(k) * dt));

            // Discard the first few steps: the filter starts with an invented
            // velocity covariance and has not yet converged, so its NIS is not
            // yet chi-square. Counting them would be measuring the initiation
            // guess, not the filter.
            if (k > 6) {
                const double d = static_cast<double>(track_nis(t, z, cfg));
                sum_nis += d;
                ++n;
                if (d < 5.991) ++under_95;
                if (d < 9.210) ++under_99;
            }
            track_update(t, z, cfg);
        }
    }

    const double mean = sum_nis / static_cast<double>(n);
    const double f95 = static_cast<double>(under_95) / static_cast<double>(n);
    const double f99 = static_cast<double>(under_99) / static_cast<double>(n);
    std::printf("       %zu samples: mean NIS %.3f (dof = 2)\n", n, mean);
    std::printf("       under the 95%% gate: %.1f%%   under the 99%%: %.1f%%\n",
                f95 * 100.0, f99 * 100.0);
    PT_CHECK_NEAR(mean, 2.0, 0.25);
    PT_CHECK_NEAR(f95, 0.95, 0.04);
    PT_CHECK_NEAR(f99, 0.99, 0.025);
}

PT_TEST(filter_beats_the_raw_measurements) {
    // The point of filtering: several noisy fixes of a predictable target beat
    // any one of them. Compared against the measurement converted straight to
    // position, which is what a tracker-less system would report.
    const TrackerConfig cfg = default_cfg();
    pt::Rng rng(5150);

    double sum_filter = 0, sum_raw = 0;
    std::size_t n = 0;
    for (std::size_t r = 0; r < 40; ++r) {
        Truth truth{100.0, 3000.0, 2.0, -5.0};
        Track t{};
        track_initiate(t, observe(truth, cfg, rng, 0), cfg, 1);
        for (std::size_t k = 1; k < 40; ++k) {
            truth.advance(1.0);
            track_predict(t, 1, cfg);
            const Measurement z = observe(truth, cfg, rng, static_cast<Real>(k));
            track_update(t, z, cfg);
            if (k < 10) continue;

            const double fx = static_cast<double>(t.state.x) - truth.x;
            const double fy = static_cast<double>(t.state.y) - truth.y;
            sum_filter += fx * fx + fy * fy;

            const double rx = static_cast<double>(z.range_m)
                            * std::sin(static_cast<double>(z.bearing_rad)) - truth.x;
            const double ry = static_cast<double>(z.range_m)
                            * std::cos(static_cast<double>(z.bearing_rad)) - truth.y;
            sum_raw += rx * rx + ry * ry;
            ++n;
        }
    }
    const double rms_filter = std::sqrt(sum_filter / static_cast<double>(n));
    const double rms_raw = std::sqrt(sum_raw / static_cast<double>(n));
    std::printf("       RMS position error: raw %.2f m, filtered %.2f m (%.2fx better)\n",
                rms_raw, rms_filter, rms_raw / rms_filter);
    PT_CHECK(rms_filter < rms_raw);
    PT_CHECK(rms_raw / rms_filter > 1.5);
}

PT_TEST(filter_recovers_velocity_a_single_detection_cannot_know) {
    // A detection gives range and bearing. Velocity comes only from watching,
    // and it is what lets a track coast through a missed scan.
    const TrackerConfig cfg = default_cfg();
    pt::Rng rng(31337);
    Truth truth{-400.0, 2000.0, 4.0, -7.0};

    Track t{};
    track_initiate(t, observe(truth, cfg, rng, 0), cfg, 1);
    for (std::size_t k = 1; k <= 40; ++k) {
        truth.advance(1.0);
        track_predict(t, 1, cfg);
        track_update(t, observe(truth, cfg, rng, static_cast<Real>(k)), cfg);
    }
    std::printf("       truth v = (%.2f, %.2f), estimate (%.2f, %.2f)\n",
                truth.vx, truth.vy, static_cast<double>(t.state.vx),
                static_cast<double>(t.state.vy));
    PT_CHECK_NEAR(t.state.vx, truth.vx, 1.5);
    PT_CHECK_NEAR(t.state.vy, truth.vy, 1.5);

    // And the closing rate, which is what the Doppler bank measures directly --
    // two independent routes to the same number.
    const double expect_rr = -(truth.x * truth.vx + truth.y * truth.vy) / truth.range();
    std::printf("       closing rate: truth %.3f m/s, tracked %.3f m/s\n",
                expect_rr, static_cast<double>(t.state.range_rate_mps()));
    PT_CHECK_NEAR(t.state.range_rate_mps(), expect_rr, 1.0);
}

// ---------------------------------------------------------------------------
// Association and track management
// ---------------------------------------------------------------------------

PT_TEST(gate_admits_the_target_and_rejects_the_distant) {
    const TrackerConfig cfg = default_cfg();
    pt::Rng rng(11);
    Truth truth{0.0, 2000.0, 0.0, 0.0};

    Track t{};
    track_initiate(t, observe(truth, cfg, rng, 0), cfg, 1);
    for (std::size_t k = 1; k < 15; ++k) {
        track_predict(t, 1, cfg);
        track_update(t, observe(truth, cfg, rng, static_cast<Real>(k)), cfg);
    }
    track_predict(t, 1, cfg);

    Measurement good;
    good.range_m = 2000;
    good.bearing_rad = 0;
    good.time_s = 15;
    PT_CHECK(static_cast<double>(track_nis(t, good, cfg))
           < static_cast<double>(cfg.gate_chi2_2dof));

    Measurement far = good;
    far.range_m = 2300;                                   // 30 sigma in range
    PT_CHECK(static_cast<double>(track_nis(t, far, cfg))
           > static_cast<double>(cfg.gate_chi2_2dof));

    Measurement off = good;
    off.bearing_rad = deg2rad(static_cast<Real>(20));     // 20 sigma in bearing
    PT_CHECK(static_cast<double>(track_nis(t, off, cfg))
           > static_cast<double>(cfg.gate_chi2_2dof));

    std::printf("       NIS: on target %.3f, 300 m off %.1f, 20 deg off %.1f (gate %.2f)\n",
                static_cast<double>(track_nis(t, good, cfg)),
                static_cast<double>(track_nis(t, far, cfg)),
                static_cast<double>(track_nis(t, off, cfg)),
                static_cast<double>(cfg.gate_chi2_2dof));

    Track dead{};
    PT_CHECK(static_cast<double>(track_nis(dead, good, cfg)) > 1e6);
    PT_CHECK(!track_update(dead, good, cfg));
}

PT_TEST(a_real_target_confirms_and_a_lone_false_alarm_does_not) {
    const TrackerConfig cfg = default_cfg();
    pt::Rng rng(909);
    clear_tracks();
    std::uint32_t next_id = 1;

    Truth truth{200.0, 1800.0, -2.0, -4.0};
    for (std::size_t k = 0; k < 8; ++k) {
        if (k > 0) truth.advance(1.0);
        std::size_t n = 0;
        g_meas[n++] = observe(truth, cfg, rng, static_cast<Real>(k));
        // One isolated false alarm on the very first scan, somewhere else.
        if (k == 0) {
            Measurement fa;
            fa.range_m = 900;
            fa.bearing_rad = deg2rad(static_cast<Real>(-25));
            fa.time_s = 0;
            g_meas[n++] = fa;
        }
        tracker_step(g_tracks, std::span<const Measurement>(g_meas.data(), n),
                     cfg, static_cast<Real>(k), next_id);
    }

    const std::size_t confirmed = count_tracks(g_tracks, TrackStatus::Confirmed);
    const std::size_t tentative = count_tracks(g_tracks, TrackStatus::Tentative);
    std::printf("       after 8 scans: %zu confirmed, %zu tentative\n",
                confirmed, tentative);
    PT_CHECK(confirmed == 1);
    PT_CHECK(tentative == 0);

    // And the surviving track is on the target.
    for (const Track& t : g_tracks) {
        if (t.status != TrackStatus::Confirmed) continue;
        PT_CHECK_NEAR(t.state.range_m(), truth.range(), 60.0);
        PT_CHECK_NEAR(rad2deg(t.state.bearing_rad()),
                      truth.bearing() * 180.0 / kPiD, 2.0);
    }
}

PT_TEST(false_alarms_do_not_survive_confirmation) {
    // The quantitative version. Poisson-ish clutter at a rate a CFAR detector
    // might produce, with no target at all: count how many confirm.
    const TrackerConfig cfg = default_cfg();
    pt::Rng rng(24601);
    clear_tracks();
    std::uint32_t next_id = 1;

    const std::size_t scans = 300;
    std::size_t total_fa = 0;
    std::size_t peak_confirmed = 0;
    for (std::size_t k = 0; k < scans; ++k) {
        std::size_t n = 0;
        // On average two false alarms per scan, uniformly scattered.
        const std::size_t count = (rng.uniform01() < 0.6) ? 2 : 1;
        for (std::size_t i = 0; i < count && n < 8; ++i) {
            Measurement fa;
            fa.range_m = static_cast<Real>(500.0 + 4500.0 * rng.uniform01());
            fa.bearing_rad = deg2rad(static_cast<Real>(-45.0 + 90.0 * rng.uniform01()));
            fa.time_s = static_cast<Real>(k);
            g_meas[n++] = fa;
        }
        total_fa += n;
        tracker_step(g_tracks, std::span<const Measurement>(g_meas.data(), n),
                     cfg, static_cast<Real>(k), next_id);
        peak_confirmed = std::max(peak_confirmed,
                                  count_tracks(g_tracks, TrackStatus::Confirmed));
    }
    std::printf("       %zu false alarms over %zu scans -> at most %zu confirmed track(s)\n",
                total_fa, scans, peak_confirmed);
    // Three-of-three confirmation on scattered clutter should essentially never
    // fire; a handful over 300 scans is the tail, not a failure.
    PT_CHECK(peak_confirmed <= 2);
}

PT_TEST(a_track_coasts_through_a_missed_scan_and_dies_after_enough) {
    const TrackerConfig cfg = default_cfg();
    pt::Rng rng(777);
    clear_tracks();
    std::uint32_t next_id = 1;

    Truth truth{0.0, 2000.0, 5.0, -3.0};
    // Six clean scans to confirm.
    for (std::size_t k = 0; k < 6; ++k) {
        if (k > 0) truth.advance(1.0);
        g_meas[0] = observe(truth, cfg, rng, static_cast<Real>(k));
        tracker_step(g_tracks, std::span<const Measurement>(g_meas.data(), 1),
                     cfg, static_cast<Real>(k), next_id);
    }
    PT_CHECK(count_tracks(g_tracks, TrackStatus::Confirmed) == 1);

    // One missed scan: the track coasts and keeps its state.
    truth.advance(1.0);
    tracker_step(g_tracks, std::span<const Measurement>(), cfg, 6, next_id);
    const std::size_t coasting = count_tracks(g_tracks, TrackStatus::Coasting);
    std::printf("       after 1 miss: %zu coasting\n", coasting);
    PT_CHECK(coasting == 1);

    // It reacquires.
    truth.advance(1.0);
    g_meas[0] = observe(truth, cfg, rng, 7);
    tracker_step(g_tracks, std::span<const Measurement>(g_meas.data(), 1), cfg, 7, next_id);
    std::printf("       after reacquisition: %zu confirmed\n",
                count_tracks(g_tracks, TrackStatus::Confirmed));
    PT_CHECK(count_tracks(g_tracks, TrackStatus::Confirmed) == 1);

    // Three misses and it is gone.
    for (std::size_t k = 8; k <= 10; ++k) {
        tracker_step(g_tracks, std::span<const Measurement>(), cfg,
                     static_cast<Real>(k), next_id);
    }
    std::printf("       after 3 misses: %zu live\n",
                g_tracks.size() - count_tracks(g_tracks, TrackStatus::Free));
    PT_CHECK(count_tracks(g_tracks, TrackStatus::Free) == g_tracks.size());
}

PT_TEST(two_targets_get_two_tracks) {
    const TrackerConfig cfg = default_cfg();
    pt::Rng rng(2024);
    clear_tracks();
    std::uint32_t next_id = 1;

    Truth a{-600.0, 2000.0, 3.0, -2.0};
    Truth b{600.0, 2600.0, -3.0, -5.0};
    for (std::size_t k = 0; k < 12; ++k) {
        if (k > 0) { a.advance(1.0); b.advance(1.0); }
        g_meas[0] = observe(a, cfg, rng, static_cast<Real>(k));
        g_meas[1] = observe(b, cfg, rng, static_cast<Real>(k));
        tracker_step(g_tracks, std::span<const Measurement>(g_meas.data(), 2),
                     cfg, static_cast<Real>(k), next_id);
    }
    const std::size_t confirmed = count_tracks(g_tracks, TrackStatus::Confirmed);
    std::printf("       two well-separated targets -> %zu confirmed tracks\n", confirmed);
    PT_CHECK(confirmed == 2);

    // Each track sits on one of them, and they have different identities.
    std::size_t matched = 0;
    for (const Track& t : g_tracks) {
        if (t.status != TrackStatus::Confirmed) continue;
        const double dr_a = std::fabs(static_cast<double>(t.state.range_m()) - a.range());
        const double dr_b = std::fabs(static_cast<double>(t.state.range_m()) - b.range());
        if (std::min(dr_a, dr_b) < 80.0) ++matched;
    }
    PT_CHECK(matched == 2);
}

// ---------------------------------------------------------------------------
// Correcting the record
// ---------------------------------------------------------------------------

PT_TEST(tracking_does_not_suppress_cross_template_ghosts) {
    // The v0.8 roadmap said time consistency would finally kill the
    // cross-template ghosts of v0.2. That was wrong, and this test says so
    // rather than the claim being quietly dropped.
    //
    // A ghost appears whenever the real arrival does, at a fixed offset set by
    // the template cross-correlation. It is therefore exactly as consistent
    // over time as the target, moves with it, and forms its own perfectly
    // healthy confirmed track. Nothing about its kinematics is wrong.
    const TrackerConfig cfg = default_cfg();
    pt::Rng rng(4711);
    clear_tracks();
    std::uint32_t next_id = 1;

    Truth truth{150.0, 2200.0, -1.0, -6.0};
    // The v0.2 measurement: an LFM arrival lights the HFM template about
    // 10.5 dB down, at a shifted lag. Here that lag offset is a fixed 120 m of
    // apparent range.
    const double ghost_offset_m = 120.0;

    for (std::size_t k = 0; k < 10; ++k) {
        if (k > 0) truth.advance(1.0);
        g_meas[0] = observe(truth, cfg, rng, static_cast<Real>(k));
        g_meas[0].label = 1;
        Truth ghost = truth;
        const double scale = (truth.range() + ghost_offset_m) / truth.range();
        ghost.x *= scale;
        ghost.y *= scale;
        ghost.vx *= scale;
        ghost.vy *= scale;
        g_meas[1] = observe(ghost, cfg, rng, static_cast<Real>(k));
        g_meas[1].label = 2;
        tracker_step(g_tracks, std::span<const Measurement>(g_meas.data(), 2),
                     cfg, static_cast<Real>(k), next_id);
    }

    const std::size_t confirmed = count_tracks(g_tracks, TrackStatus::Confirmed);
    std::printf("       one target + its cross-template ghost -> %zu confirmed tracks\n",
                confirmed);
    std::printf("       The ghost tracks as cleanly as the target: it arrives whenever\n");
    std::printf("       the target does, at a fixed offset, so time consistency has\n");
    std::printf("       nothing to object to. Tracking kills FALSE ALARMS, which do not\n");
    std::printf("       repeat -- not ghosts, which do. Suppressing these needs the\n");
    std::printf("       offset and amplitude ratio to be recognised as a template\n");
    std::printf("       artefact, which this library does not do.\n");
    PT_CHECK(confirmed == 2);   // the honest outcome, asserted so it stays true
}

PT_TEST(tracker_handles_degenerate_input) {
    const TrackerConfig cfg = default_cfg();
    clear_tracks();
    std::uint32_t next_id = 1;

    // No measurements, no tracks: nothing happens and nothing breaks.
    PT_CHECK(tracker_step(g_tracks, std::span<const Measurement>(), cfg, 0, next_id) == 0);
    PT_CHECK(tracker_step(std::span<Track>(),
                          std::span<const Measurement>(g_meas.data(), 0), cfg, 0, next_id) == 0);

    // A measurement at the origin has no bearing; the track initiates but the
    // filter refuses to update it rather than dividing by zero.
    Measurement z;
    z.range_m = 0;
    z.bearing_rad = 0;
    z.time_s = 0;
    Track t{};
    track_initiate(t, z, cfg, 1);
    PT_CHECK(t.status == TrackStatus::Tentative);
    PT_CHECK(!track_update(t, z, cfg));
    PT_CHECK(static_cast<double>(track_nis(t, z, cfg)) > 1e6);

    // A dead track is not predicted.
    Track dead{};
    const Real before = dead.state.x;
    track_predict(dead, 1, cfg);
    PT_CHECK(dead.state.x == before);
    PT_CHECK(dead.age == 0);
}

// ---------------------------------------------------------------------------
// v0.10: fusing the range rate the Doppler bank already measures
// ---------------------------------------------------------------------------

PT_TEST(range_rate_fusion_keeps_the_filter_consistent) {
    // Adding a third measurement changes the NIS distribution from 2 dof to 3.
    // If the new Jacobian row or the new R entry is wrong the filter still
    // tracks -- and this is the only thing that notices.
    TrackerConfig cfg = default_cfg();
    cfg.range_rate_sigma_mps = 2;
    cfg.gate_chi2_3dof = chi2_gate(static_cast<Real>(0.99), 3);
    pt::Rng rng(880808);

    double sum_nis = 0;
    std::size_t n = 0, under_95 = 0, under_99 = 0;
    const double g95 = static_cast<double>(chi2_gate(static_cast<Real>(0.95), 3));
    const double g99 = static_cast<double>(chi2_gate(static_cast<Real>(0.99), 3));

    for (std::size_t r = 0; r < 60; ++r) {
        Truth truth{-150.0 + 10.0 * rng.uniform(), 2400.0 + 300.0 * rng.uniform(),
                     3.0 * rng.uniform(), -6.0 + 2.0 * rng.uniform()};
        Track t{};
        Measurement z0 = observe(truth, cfg, rng, 0);
        track_initiate(t, z0, cfg, 1);

        for (std::size_t k = 1; k < 40; ++k) {
            truth.advance(1.0);
            track_predict(t, 1, cfg);
            Measurement z = observe(truth, cfg, rng, static_cast<Real>(k));
            const double true_rr = -(truth.x * truth.vx + truth.y * truth.vy) / truth.range();
            z.range_rate_mps = static_cast<Real>(
                true_rr + static_cast<double>(cfg.range_rate_sigma_mps) * rng.normal());
            z.has_range_rate = true;

            if (k > 6) {
                const double d = static_cast<double>(track_nis(t, z, cfg));
                sum_nis += d;
                ++n;
                if (d < g95) ++under_95;
                if (d < g99) ++under_99;
            }
            track_update(t, z, cfg);
        }
    }
    const double mean = sum_nis / static_cast<double>(n);
    std::printf("       3-dof gates: 95%% at %.3f, 99%% at %.3f\n", g95, g99);
    std::printf("       %zu samples: mean NIS %.3f (dof = 3)\n", n, mean);
    std::printf("       under the 95%% gate: %.1f%%   under the 99%%: %.1f%%\n",
                100.0 * static_cast<double>(under_95) / static_cast<double>(n),
                100.0 * static_cast<double>(under_99) / static_cast<double>(n));
    PT_CHECK_NEAR(g95, 7.815, 0.01);
    PT_CHECK_NEAR(g99, 11.345, 0.02);
    PT_CHECK_NEAR(mean, 3.0, 0.35);
    PT_CHECK_NEAR(static_cast<double>(under_95) / static_cast<double>(n), 0.95, 0.04);
}

PT_TEST(range_rate_fusion_sharpens_the_velocity_estimate) {
    // The payoff. Position history infers velocity slowly; a direct measurement
    // of its radial component observes it immediately. The improvement should
    // be largest early, when the position history is shortest.
    TrackerConfig cfg = default_cfg();
    cfg.range_rate_sigma_mps = 2;

    auto run = [&](bool fuse, std::size_t steps) {
        pt::Rng rng(1357);
        double sum_sq = 0;
        const std::size_t runs = 60;
        for (std::size_t r = 0; r < runs; ++r) {
            Truth truth{200.0, 2500.0, 2.0, -6.0};
            Track t{};
            track_initiate(t, observe(truth, cfg, rng, 0), cfg, 1);
            for (std::size_t k = 1; k <= steps; ++k) {
                truth.advance(1.0);
                track_predict(t, 1, cfg);
                Measurement z = observe(truth, cfg, rng, static_cast<Real>(k));
                if (fuse) {
                    const double rr = -(truth.x * truth.vx + truth.y * truth.vy) / truth.range();
                    z.range_rate_mps = static_cast<Real>(
                        rr + static_cast<double>(cfg.range_rate_sigma_mps) * rng.normal());
                    z.has_range_rate = true;
                }
                track_update(t, z, cfg);
            }
            // RADIAL velocity error. A range rate observes only the component
            // along the line of sight, so that is the only component it can
            // improve -- measuring total velocity error would dilute the
            // effect with a tangential term the measurement never sees.
            const double true_rr = -(truth.x * truth.vx + truth.y * truth.vy) / truth.range();
            const double e = static_cast<double>(t.state.range_rate_mps()) - true_rr;
            sum_sq += e * e;
        }
        return std::sqrt(sum_sq / static_cast<double>(runs));
    };

    std::printf("       radial velocity error (the only component a range rate observes)\n");
    std::printf("       %8s %16s %16s %10s\n",
                "scans", "position only", "with range rate", "better by");
    for (std::size_t steps : {std::size_t(3), std::size_t(6), std::size_t(20)}) {
        const double a = run(false, steps);
        const double b = run(true, steps);
        std::printf("       %8zu %14.3f m/s %14.3f m/s %9.2fx\n", steps, a, b, a / b);
        // Fusion must never make things materially worse, but by 20 scans the
        // filter's own estimate (0.37 m/s) is already better than the 2 m/s
        // measurement, so there is nothing left for it to add. A measurement
        // helps exactly while it beats the estimate you already have.
        PT_CHECK(b < a * 1.05);
    }
    // The gain is largest early, when the position history is shortest and the
    // filter has least to infer velocity from.
    const double early = run(false, 3) / run(true, 3);
    const double late = run(false, 20) / run(true, 20);
    std::printf("       gain at 3 scans %.2fx, at 20 scans %.2fx -- it fades as the\n"
                "       position history grows long enough to infer velocity itself\n",
                early, late);
    PT_CHECK(early > 1.5);
    PT_CHECK(early > late);
}

PT_TEST(an_unresolved_doppler_bin_is_not_a_measurement_of_zero) {
    // The PulseDescriptor reports 0 m/s when the bank has no Doppler coverage.
    // Feeding that in as a measurement would pin every track's radial velocity
    // to zero, so `has_range_rate` exists to say "the bank did not resolve it".
    TrackerConfig cfg = default_cfg();
    pt::Rng rng(246);
    Truth truth{0.0, 2000.0, 0.0, -8.0};   // closing hard at 8 m/s

    Track fused{}, unfused{};
    track_initiate(fused, observe(truth, cfg, rng, 0), cfg, 1);
    track_initiate(unfused, observe(truth, cfg, rng, 0), cfg, 2);

    for (std::size_t k = 1; k <= 15; ++k) {
        truth.advance(1.0);
        track_predict(fused, 1, cfg);
        track_predict(unfused, 1, cfg);
        const Measurement base = observe(truth, cfg, rng, static_cast<Real>(k));

        // The wrong thing: a zero from an unresolved bin, taken at face value.
        Measurement wrong = base;
        wrong.range_rate_mps = 0;
        wrong.has_range_rate = true;
        track_update(fused, wrong, cfg);

        // The right thing: the flag left clear, so the filter ignores it.
        track_update(unfused, base, cfg);
    }
    std::printf("       truth closing 8.00 m/s\n");
    std::printf("       zero taken as a measurement -> %.3f m/s\n",
                static_cast<double>(fused.state.range_rate_mps()));
    std::printf("       flag left clear             -> %.3f m/s\n",
                static_cast<double>(unfused.state.range_rate_mps()));
    // Believing the zero drags the estimate towards it; ignoring it does not.
    PT_CHECK(static_cast<double>(fused.state.range_rate_mps())
           < static_cast<double>(unfused.state.range_rate_mps()));
    PT_CHECK_NEAR(unfused.state.range_rate_mps(), 8.0, 1.5);
}

// ---------------------------------------------------------------------------
// v0.10: ghost recognition
// ---------------------------------------------------------------------------

namespace {

// Runs a two-object scenario and returns how many confirmed tracks survive
// ghost suppression. `same_label` decides whether the second object returns the
// same waveform as the first -- which is what separates a real formation from a
// template artefact.
std::size_t run_pair(double offset_m, double bearing_offset_deg, bool same_label,
                     double amplitude_ratio, std::uint64_t seed) {
    const TrackerConfig cfg = default_cfg();
    GhostConfig ghost;
    // Set the pairing tolerance from the bearing accuracy the array actually
    // delivers -- three sigma of the measurement. Two tracks on one true
    // bearing routinely differ by more than one sigma, so a tolerance chosen
    // from the geometry rather than the noise never fires.
    ghost.max_bearing_delta_rad = static_cast<Real>(3) * cfg.bearing_sigma_rad;
    pt::Rng rng(seed);
    clear_tracks();
    std::uint32_t next_id = 1;

    Truth truth{150.0, 2200.0, -1.0, -6.0};
    for (std::size_t k = 0; k < 20; ++k) {
        if (k > 0) truth.advance(1.0);
        g_meas[0] = observe(truth, cfg, rng, static_cast<Real>(k));
        g_meas[0].label = 1;
        g_meas[0].amplitude = 1;

        Truth second = truth;
        const double scale = (truth.range() + offset_m) / truth.range();
        second.x *= scale;
        second.y *= scale;
        second.vx *= scale;
        second.vy *= scale;
        g_meas[1] = observe(second, cfg, rng, static_cast<Real>(k));
        g_meas[1].bearing_rad += deg2rad(static_cast<Real>(bearing_offset_deg));
        g_meas[1].label = same_label ? 1 : 2;
        g_meas[1].amplitude = static_cast<Real>(amplitude_ratio);

        tracker_step(g_tracks, std::span<const Measurement>(g_meas.data(), 2),
                     cfg, static_cast<Real>(k), next_id);
    }
    suppress_template_ghosts(g_tracks, ghost);
    // Confirmed OR Coasting: a track that happened to miss the last scan is
    // still a track, and counting only Confirmed would report a detector's
    // miss rate rather than the number of contacts.
    return count_established(g_tracks);
}

}  // namespace

PT_TEST(ghost_recognition_removes_the_artefact) {
    // v0.9 measured that tracking alone leaves two tracks here. Recognising the
    // ghost as a template artefact -- same bearing, same kinematics, weaker,
    // DIFFERENT waveform -- removes it.
    const std::size_t n = run_pair(120.0, 0.0, /*same_label=*/false, 0.30, 4711);
    std::printf("       target + cross-template ghost -> %zu established track(s)\n", n);
    PT_CHECK(n == 1);

    // And the survivor is the strong one, at the target's range.
    for (const Track& t : g_tracks) {
        if (t.status != TrackStatus::Confirmed && t.status != TrackStatus::Coasting) continue;
        PT_CHECK(t.amplitude > static_cast<Real>(0.8));
        PT_CHECK(t.label == 1);
    }
}

PT_TEST(ghost_recognition_spares_genuine_targets) {
    // The discriminating direction. Suppression must not eat real contacts.
    // A second target at a different bearing: kept.
    PT_CHECK(run_pair(120.0, 6.0, false, 0.30, 991) == 2);
    // A second target of comparable strength: kept, whatever its bearing.
    PT_CHECK(run_pair(120.0, 0.0, false, 0.95, 992) == 2);
    // A second target far beyond any template cross-correlation lag: kept.
    PT_CHECK(run_pair(1500.0, 0.0, false, 0.30, 993) == 2);
    std::printf("       different bearing / comparable strength / distant: all kept\n");
}

PT_TEST(the_label_check_is_what_saves_a_line_astern_formation) {
    // The honest failure this design avoids, demonstrated both ways.
    //
    // Two real targets in line astern -- same bearing, same course and speed,
    // fixed separation, one weaker because it is further away -- are
    // kinematically INDISTINGUISHABLE from a ghost pair. The only thing that
    // separates them is the waveform: a ghost is by definition a different
    // template, while two targets lit by one sonar return the same one.
    const std::size_t with_label_check = run_pair(120.0, 0.0, /*same_label=*/true, 0.30, 5150);
    std::printf("       line-astern formation, same waveform -> %zu tracks (kept)\n",
                with_label_check);
    PT_CHECK(with_label_check == 2);

    // The same geometry with different waveforms is suppressed -- which is the
    // right call for a ghost and the wrong one for two targets that happen to
    // return different types. That residual failure mode is real and is
    // documented rather than papered over.
    const std::size_t without = run_pair(120.0, 0.0, false, 0.30, 5150);
    std::printf("       same geometry, different waveforms    -> %zu track (suppressed)\n",
                without);
    PT_CHECK(without == 1);
}

// ---------------------------------------------------------------------------
// v0.10: association in global cost order
// ---------------------------------------------------------------------------

PT_TEST(crossing_targets_do_not_swap_tracks) {
    // Per-track greedy assigns in TRACK order, so whichever track is examined
    // first takes the measurement it likes best -- even when the other track
    // wants it far more. During a crossing that is exactly how identities swap.
    const TrackerConfig cfg = default_cfg();
    pt::Rng rng(31415);
    clear_tracks();
    std::uint32_t next_id = 1;

    // Two targets converging on the same point from opposite bearings.
    Truth a{-400.0, 2000.0, 20.0, 0.0};
    Truth b{ 400.0, 2000.0, -20.0, 0.0};

    std::uint32_t id_a = 0, id_b = 0;
    for (std::size_t k = 0; k < 40; ++k) {
        if (k > 0) { a.advance(1.0); b.advance(1.0); }
        g_meas[0] = observe(a, cfg, rng, static_cast<Real>(k));
        g_meas[1] = observe(b, cfg, rng, static_cast<Real>(k));
        tracker_step(g_tracks, std::span<const Measurement>(g_meas.data(), 2),
                     cfg, static_cast<Real>(k), next_id);

        // Latch the identities once both are confirmed and still well apart.
        if (k == 8) {
            for (const Track& t : g_tracks) {
                if (t.status != TrackStatus::Confirmed) continue;
                if (t.state.x < 0) id_a = t.id; else id_b = t.id;
            }
        }
    }
    PT_CHECK(id_a != 0);
    PT_CHECK(id_b != 0);

    // After the crossing, the track that started on the left must be the one
    // now on the right, and keep its identity through it.
    std::uint32_t now_left = 0, now_right = 0;
    for (const Track& t : g_tracks) {
        if (t.status != TrackStatus::Confirmed) continue;
        if (t.state.x < 0) now_left = t.id; else now_right = t.id;
    }
    std::printf("       before crossing: left=%u right=%u\n", id_a, id_b);
    std::printf("       after  crossing: left=%u right=%u\n", now_left, now_right);
    PT_CHECK(count_tracks(g_tracks, TrackStatus::Confirmed) == 2);
    // a started left moving right, so it ends on the right.
    PT_CHECK(now_right == id_a);
    PT_CHECK(now_left == id_b);
}
