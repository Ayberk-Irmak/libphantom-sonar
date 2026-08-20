// SPDX-License-Identifier: Apache-2.0
// Coordinated turn with the turn rate estimated, and the IMM tracker.
//
// Two claims to test. First, that estimating omega actually buys what
// bracketing it cannot: a rate outside the bracket. Second, that it costs
// something, because a fifth state always does -- and where.
#include "framework.hpp"

#include "phantom/ctrv.hpp"
#include "phantom/imm.hpp"
#include "phantom/tracker.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>

using namespace phantom;

namespace {

constexpr Real kDt = 1;

struct Truth { Real x, y, vx, vy; };

void advance(Truth& t, Real omega, Real dt) {
    if (std::fabs(static_cast<double>(omega)) < 1e-12) {
        t.x += t.vx * dt; t.y += t.vy * dt;
        return;
    }
    const Real c = std::cos(omega * dt), s = std::sin(omega * dt);
    const Real nx = t.x + (s / omega) * t.vx - ((1 - c) / omega) * t.vy;
    const Real ny = t.y + ((1 - c) / omega) * t.vx + (s / omega) * t.vy;
    const Real nvx = c * t.vx - s * t.vy;
    const Real nvy = s * t.vx + c * t.vy;
    t.x = nx; t.y = ny; t.vx = nvx; t.vy = nvy;
}

Measurement measure(const Truth& t, Real time_s, pt::Rng& rng, const TrackerConfig& cfg) {
    Measurement z;
    const double r = std::sqrt(static_cast<double>(t.x * t.x + t.y * t.y));
    const double b = std::atan2(static_cast<double>(t.x), static_cast<double>(t.y));
    z.range_m = static_cast<Real>(r + static_cast<double>(cfg.range_sigma_m) * rng.normal());
    z.bearing_rad = static_cast<Real>(b + static_cast<double>(cfg.bearing_sigma_rad) * rng.normal());
    z.time_s = time_s;
    return z;
}

double pos_err(Real sx, Real sy, const Truth& t) {
    const double dx = static_cast<double>(sx - t.x), dy = static_cast<double>(sy - t.y);
    return std::sqrt(dx * dx + dy * dy);
}

TrackerConfig sonar_config() {
    TrackerConfig cfg;
    cfg.range_sigma_m = 5;
    cfg.bearing_sigma_rad = deg2rad(static_cast<Real>(0.5));
    cfg.process_accel_mps2 = static_cast<Real>(1.0);
    return cfg;
}

}  // namespace

PT_TEST(ctrv_recovers_a_turn_rate_the_imm_can_only_saturate_at) {
    // The reason a fifth state is worth its cost. Both filters are given
    // manoeuvre models set for 3 deg/s; the target turns at 8. The IMM's
    // estimate cannot leave [-3, +3] by construction, so it saturates and
    // reports a turn three times gentler than the one happening.
    const TrackerConfig cfg = sonar_config();
    ImmConfig imm;
    imm.turn_rate_rps = deg2rad(static_cast<Real>(3));
    CtrvConfig ct;

    std::printf("       models bracket at 3 deg/s; measuring what each reports\n");
    std::printf("       %10s %14s %16s\n", "truth", "IMM (bracketed)", "CTRV (estimated)");
    for (int truth_deg : {2, 5, 8}) {
        pt::Rng rng(5150 + static_cast<std::uint64_t>(truth_deg));
        Truth truth{2000, 9000, -35, 5};
        const Real omega = deg2rad(static_cast<Real>(truth_deg));

        ImmTrack it;
        CtrvTrack ctk;
        const Measurement z0 = measure(truth, 0, rng, cfg);
        imm_initiate(it, z0, cfg, imm, 1);
        ctrv_initiate(ctk, z0, cfg, ct, 1);

        for (std::size_t k = 1; k <= 40; ++k) {
            advance(truth, omega, kDt);
            const Measurement z = measure(truth, static_cast<Real>(k) * kDt, rng, cfg);
            imm_predict(it, kDt, cfg, imm);
            imm_update(it, z, cfg, imm);
            ctrv_predict(ctk, kDt, cfg, ct);
            ctrv_update(ctk, z, cfg, ct);
        }
        const double imm_w = static_cast<double>(rad2deg(imm_turn_rate_estimate(it, imm)));
        const double ctrv_w = static_cast<double>(rad2deg(ctk.state.turn_rate_rps));
        std::printf("       %8d d/s %11.2f d/s %13.2f d/s\n", truth_deg, imm_w, ctrv_w);

        // CTRV must be the closer of the two whenever the truth is outside the
        // bracket, and must not saturate.
        if (truth_deg > 3) {
            PT_CHECK(std::fabs(ctrv_w - truth_deg) < std::fabs(imm_w - truth_deg));
            PT_CHECK(ctrv_w > 3.0);
        }
        PT_CHECK(std::fabs(imm_w) <= 3.0 + 1e-6);   // the bracket, by construction
    }
}

PT_TEST(turn_rate_process_noise_is_the_knob_and_zero_is_the_trap) {
    // The fifth state's process noise decides whether omega can be LEARNED, and
    // setting it to zero fails in the way that is hardest to notice: the filter
    // converges, reports a small uncertainty, and is wrong.
    //
    // With no random walk on omega the covariance shrinks monotonically, the
    // Kalman gain on the fifth state goes to nothing, and the filter stops being
    // able to follow a change. It then reports 0.1 deg/s of uncertainty about a
    // number that is out by 3.
    const TrackerConfig cfg = sonar_config();
    std::printf("       target flies straight for 30 scans, then turns at 5.00 deg/s\n");
    std::printf("       %10s %16s %14s\n", "q_w", "reported sigma", "estimate");

    double sigma_at[4] = {0, 0, 0, 0};
    double est_at[4] = {0, 0, 0, 0};
    const double qs[4] = {0.0, 0.005, 0.01, 0.02};
    for (int i = 0; i < 4; ++i) {
        CtrvConfig ct;
        ct.turn_rate_noise_rps = static_cast<Real>(qs[i]);
        pt::Rng rng(8080);
        Truth truth{1500, 7000, -25, 8};
        CtrvTrack t;
        ctrv_initiate(t, measure(truth, 0, rng, cfg), cfg, ct, 1);
        for (std::size_t k = 1; k <= 60; ++k) {
            const Real w = (k > 30) ? deg2rad(static_cast<Real>(5)) : static_cast<Real>(0);
            advance(truth, w, kDt);
            const Measurement z = measure(truth, static_cast<Real>(k) * kDt, rng, cfg);
            ctrv_predict(t, kDt, cfg, ct);
            ctrv_update(t, z, cfg, ct);
        }
        sigma_at[i] = static_cast<double>(rad2deg(ctrv_turn_rate_sigma_rps(t)));
        est_at[i] = static_cast<double>(rad2deg(t.state.turn_rate_rps));
        std::printf("       %10.3f %13.2f deg/s %10.2f deg/s%s\n",
                    qs[i], sigma_at[i], est_at[i],
                    (i == 0) ? "   <- converged, and wrong" : "");
    }

    // q_w = 0 is confidently wrong: the smallest reported uncertainty of the
    // four, and the largest error.
    PT_CHECK(sigma_at[0] < sigma_at[1]);
    PT_CHECK(std::fabs(est_at[0] - 5.0) > 2.0);
    // Every non-zero setting recovers the rate.
    for (int i = 1; i < 4; ++i) PT_CHECK(std::fabs(est_at[i] - 5.0) < 1.0);
    // And more process noise buys a larger reported uncertainty, monotonically.
    PT_CHECK(sigma_at[1] < sigma_at[2]);
    PT_CHECK(sigma_at[2] < sigma_at[3]);

    std::printf("       The reported sigma reaches a steady state rather than collapsing,\n");
    std::printf("       which is correct: omega can change at any moment, so a filter that\n");
    std::printf("       stops allowing for that has stopped being a tracker.\n");
}

PT_TEST(ctrv_pays_for_the_fifth_state_on_a_straight_target) {
    // The cost, measured rather than waved at. An extra state is an extra thing
    // to get wrong, and a straight target gives it no information -- so the
    // filter spends accuracy on estimating something that is not happening.
    const TrackerConfig cfg = sonar_config();
    CtrvConfig ct;
    double ctrv_sse = 0, ekf_sse = 0;
    std::size_t n = 0;
    for (std::uint64_t seed = 0; seed < 12; ++seed) {
        pt::Rng rng(9000 + seed);
        Truth truth{2000, 8000, -30, 10};
        CtrvTrack c;
        Track e;
        const Measurement z0 = measure(truth, 0, rng, cfg);
        ctrv_initiate(c, z0, cfg, ct, 1);
        track_initiate(e, z0, cfg, 1);
        for (std::size_t k = 1; k <= 50; ++k) {
            advance(truth, 0, kDt);
            const Measurement z = measure(truth, static_cast<Real>(k) * kDt, rng, cfg);
            ctrv_predict(c, kDt, cfg, ct);
            ctrv_update(c, z, cfg, ct);
            track_predict(e, kDt, cfg);
            track_update(e, z, cfg);
            if (k < 6) continue;
            const double a = pos_err(c.state.x, c.state.y, truth);
            const double b = pos_err(e.state.x, e.state.y, truth);
            ctrv_sse += a * a; ekf_sse += b * b; ++n;
        }
    }
    const double ctrv_rms = std::sqrt(ctrv_sse / static_cast<double>(n));
    const double ekf_rms = std::sqrt(ekf_sse / static_cast<double>(n));
    std::printf("       straight target, 12 runs: CTRV %.2f m, plain CV EKF %.2f m (%.2fx)\n",
                ctrv_rms, ekf_rms, ctrv_rms / ekf_rms);
    std::printf("       The fifth state is not free. Use CTRV when you need the RATE;\n");
    std::printf("       use the IMM when you need position through a manoeuvre.\n");
    // It must not be catastrophic -- a filter that falls apart on the easy case
    // is not usable -- but it is allowed to be worse, and it is.
    PT_CHECK(ctrv_rms < ekf_rms * 1.6);
}

PT_TEST(the_ctrv_zero_turn_limit_matches_constant_velocity) {
    // Same 0/0 as the IMM's transition, but here it appears in the JACOBIAN too,
    // and the series threshold must be the same in both -- a state and a
    // covariance that switch at different points describe different filters.
    const TrackerConfig cfg = sonar_config();
    double previous = 1e30;
    for (int e = 2; e <= 8; ++e) {
        CtrvConfig ct;
        ct.turn_rate_noise_rps = 0;          // isolate the transition
        ct.process_accel_mps2 = 0;
        Measurement z;
        z.range_m = 6000;
        z.bearing_rad = deg2rad(static_cast<Real>(15));

        CtrvTrack turning, straight;
        ctrv_initiate(turning, z, cfg, ct, 1);
        ctrv_initiate(straight, z, cfg, ct, 2);
        turning.state.vx = 14;  turning.state.vy = -9;
        straight.state.vx = 14; straight.state.vy = -9;
        turning.state.turn_rate_rps = static_cast<Real>(std::pow(10.0, -e));
        straight.state.turn_rate_rps = 0;

        ctrv_predict(turning, 10, cfg, ct);
        ctrv_predict(straight, 10, cfg, ct);
        const double d = pos_err(turning.state.x, turning.state.y,
                                 Truth{straight.state.x, straight.state.y, 0, 0});
        std::printf("       omega = 1e-%d rad/s: |CT - CV| after 10 s = %.3e m\n", e, d);
        const double floor_m = pt::tol(1e-9, 1e-3);
        if (previous > floor_m) PT_CHECK(d < previous);
        previous = d;
    }
    PT_CHECK(previous < pt::tol(1e-5, 1e-3));
}

PT_TEST(imm_tracker_step_manages_tracks_like_the_single_model_one) {
    // The integration v0.11 deliberately left out. Same association, same
    // M-of-N management, gating on the combined estimate.
    const TrackerConfig cfg = sonar_config();
    ImmConfig imm;
    pt::Rng rng(1234);

    std::array<ImmTrack, 16> tracks{};
    std::uint32_t next_id = 1;
    Truth a{-1500, 9000, 25, -6};
    Truth b{ 1500, 9000, -25, -6};

    std::size_t confirmed_at_end = 0;
    for (std::size_t k = 1; k <= 40; ++k) {
        // Both targets turn, in opposite directions, from scan 15.
        advance(a, (k >= 15) ? deg2rad(static_cast<Real>(4)) : static_cast<Real>(0), kDt);
        advance(b, (k >= 15) ? deg2rad(static_cast<Real>(-4)) : static_cast<Real>(0), kDt);
        std::array<Measurement, 2> zs{measure(a, static_cast<Real>(k), rng, cfg),
                                      measure(b, static_cast<Real>(k), rng, cfg)};
        imm_tracker_step(tracks, zs, cfg, imm, static_cast<Real>(k), next_id);
        if (k == 40) confirmed_at_end = imm_count_established(tracks);
    }
    std::printf("       two turning targets, 40 scans -> %zu established track(s)\n",
                confirmed_at_end);
    PT_CHECK(confirmed_at_end == 2);

    // Each track should have noticed its own manoeuvre, in the right direction.
    int signs = 0;
    for (const ImmTrack& t : tracks) {
        if (t.status != TrackStatus::Confirmed && t.status != TrackStatus::Coasting) continue;
        const double w = static_cast<double>(rad2deg(imm_turn_rate_estimate(t, imm)));
        std::printf("       track %u: P(manoeuvre) %.3f, turn %+.2f deg/s\n",
                    t.id, static_cast<double>(imm_manoeuvre_probability(t)), w);
        PT_CHECK(imm_manoeuvre_probability(t) > static_cast<Real>(0.5));
        signs += (w > 0) ? 1 : -1;
    }
    // One turned each way, so the signs must cancel.
    PT_CHECK(signs == 0);
}

PT_TEST(imm_tracker_step_rejects_false_alarms_like_the_single_model_one) {
    // M-of-N confirmation must survive the swap to a mixture. Pure clutter,
    // no target at all.
    const TrackerConfig cfg = sonar_config();
    ImmConfig imm;
    pt::Rng rng(60006);
    std::array<ImmTrack, 32> tracks{};
    std::uint32_t next_id = 1;
    std::size_t total_false = 0, ever_confirmed = 0;

    for (std::size_t k = 1; k <= 200; ++k) {
        std::array<Measurement, 3> zs{};
        std::size_t n = 0;
        for (std::size_t j = 0; j < 3; ++j) {
            if (rng.uniform01() > 0.6) continue;
            Measurement z;
            z.range_m = static_cast<Real>(3000.0 + 6000.0 * rng.uniform01());
            z.bearing_rad = deg2rad(static_cast<Real>(-30.0 + 60.0 * rng.uniform01()));
            z.time_s = static_cast<Real>(k);
            zs[n++] = z;
            ++total_false;
        }
        imm_tracker_step(tracks, std::span<const Measurement>(zs.data(), n),
                         cfg, imm, static_cast<Real>(k), next_id);
        ever_confirmed = std::max(ever_confirmed, imm_count_established(tracks));
    }
    std::printf("       %zu false alarms over 200 scans -> at most %zu established track(s)\n",
                total_false, ever_confirmed);
    PT_CHECK(ever_confirmed == 0);
}
