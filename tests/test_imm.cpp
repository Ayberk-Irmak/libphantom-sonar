// SPDX-License-Identifier: Apache-2.0
// Interacting multiple model filter.
//
// The claim to verify is not "the IMM works" but the specific trade it buys:
// it should beat a single CV filter on a MANOEUVRING target without paying for
// it on a straight one. A filter that only wins on turns is a filter that has
// been given more process noise, and that is not what an IMM is for.
#include "framework.hpp"

#include "phantom/imm.hpp"
#include "phantom/tracker.hpp"

#include <cmath>
#include <cstdio>

using namespace phantom;

namespace {

constexpr Real kDt = 1;

struct Truth { Real x, y, vx, vy; };

// Advances truth by a coordinated turn of `omega` rad/s.
void advance(Truth& t, Real omega, Real dt) {
    if (std::fabs(static_cast<double>(omega)) < 1e-12) {
        t.x += t.vx * dt;
        t.y += t.vy * dt;
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

double position_error(const TargetState& s, const Truth& t) {
    const double dx = static_cast<double>(s.x - t.x);
    const double dy = static_cast<double>(s.y - t.y);
    return std::sqrt(dx * dx + dy * dy);
}

// Runs one scenario through both filters. `turn_start`/`turn_end` bound a
// coordinated turn at `omega`; outside them the target flies straight.
struct Outcome {
    double imm_rms, ekf_rms, peak_manoeuvre_prob;
    // Split by phase, because the aggregate hides where the gain comes from.
    double imm_during, ekf_during;   // while the target is turning
    double imm_after,  ekf_after;    // after it straightens out again
};

Outcome run(Real omega, std::size_t turn_start, std::size_t turn_end,
            std::size_t scans, std::uint64_t seed, bool trace) {
    TrackerConfig cfg;
    cfg.range_sigma_m = 5;
    cfg.bearing_sigma_rad = deg2rad(static_cast<Real>(0.5));
    // The single-model filter gets process noise BETWEEN the IMM's two, so it
    // is not being handicapped: it is tuned as well as one number allows.
    cfg.process_accel_mps2 = static_cast<Real>(1.0);
    ImmConfig imm;
    imm.turn_rate_rps = deg2rad(static_cast<Real>(3));

    pt::Rng rng(seed);
    Truth truth{2000, 8000, -30, 10};

    ImmTrack it;
    Track ek;
    Measurement z0 = measure(truth, 0, rng, cfg);
    imm_initiate(it, z0, cfg, imm);
    track_initiate(ek, z0, cfg, 1);

    double imm_sse = 0, ekf_sse = 0, peak = 0;
    double imm_d = 0, ekf_d = 0, imm_a = 0, ekf_a = 0;
    std::size_t counted = 0, n_d = 0, n_a = 0;
    for (std::size_t k = 1; k <= scans; ++k) {
        const Real w = (k >= turn_start && k < turn_end) ? omega : static_cast<Real>(0);
        advance(truth, w, kDt);
        const Measurement z = measure(truth, static_cast<Real>(k) * kDt, rng, cfg);

        imm_predict(it, kDt, cfg, imm);
        imm_update(it, z, cfg, imm);
        track_predict(ek, kDt, cfg);
        track_update(ek, z, cfg);

        const double pm = static_cast<double>(imm_manoeuvre_probability(it));
        peak = std::max(peak, pm);
        if (trace && (k % 4 == 0 || k == turn_start || k == turn_end)) {
            std::printf("       k=%2zu %s  P(manoeuvre)=%.3f  w_est=%+6.2f deg/s"
                        "  IMM %6.1f m  EKF %6.1f m\n",
                        k, (w != 0 ? "TURN " : "     "), pm,
                        static_cast<double>(rad2deg(imm_turn_rate_estimate(it, imm))),
                        position_error(it.state, truth), position_error(ek.state, truth));
        }
        // Skip the first few scans: both filters are still converging from a
        // single detection and the comparison says nothing about the models.
        if (k < 5) continue;
        const double ei = position_error(it.state, truth);
        const double ee = position_error(ek.state, truth);
        imm_sse += ei * ei;
        ekf_sse += ee * ee;
        ++counted;
        if (k >= turn_start && k < turn_end) { imm_d += ei * ei; ekf_d += ee * ee; ++n_d; }
        else if (k >= turn_end)              { imm_a += ei * ei; ekf_a += ee * ee; ++n_a; }
    }
    auto rms = [](double sse, std::size_t n) {
        return n ? std::sqrt(sse / static_cast<double>(n)) : 0.0;
    };
    return Outcome{rms(imm_sse, counted), rms(ekf_sse, counted), peak,
                   rms(imm_d, n_d), rms(ekf_d, n_d), rms(imm_a, n_a), rms(ekf_a, n_a)};
}

}  // namespace

PT_TEST(imm_beats_a_single_model_on_a_manoeuvring_target) {
    std::printf("       target turns at 3 deg/s between scans 15 and 35:\n");
    const Outcome o = run(deg2rad(static_cast<Real>(3)), 15, 35, 60, 4242, true);
    std::printf("       RMS position error       IMM      EKF    ratio\n");
    std::printf("         during the turn     %6.2f m %6.2f m   %.2fx\n",
                o.imm_during, o.ekf_during, o.ekf_during / o.imm_during);
    std::printf("         after it ends       %6.2f m %6.2f m   %.2fx\n",
                o.imm_after, o.ekf_after, o.ekf_after / o.imm_after);
    std::printf("         overall             %6.2f m %6.2f m   %.2fx\n",
                o.imm_rms, o.ekf_rms, o.ekf_rms / o.imm_rms);
    std::printf("       The gain is mostly in RECOVERY, not in the turn. During the\n");
    std::printf("       manoeuvre both filters lag and the single model's larger process\n");
    std::printf("       noise partly covers for it; afterwards the IMM switches back to a\n");
    std::printf("       quiet model within a scan or two while the EKF is still coasting\n");
    std::printf("       on the velocity the turn left it with. That asymmetry is the\n");
    std::printf("       honest description of what an IMM buys here.\n");
    PT_CHECK(o.imm_rms < o.ekf_rms);
    PT_CHECK(o.imm_after < o.ekf_after);
    // And the manoeuvre is DETECTED, not merely survived.
    std::printf("       peak P(manoeuvre) during the turn: %.3f\n", o.peak_manoeuvre_prob);
    PT_CHECK(o.peak_manoeuvre_prob > 0.8);
}

PT_TEST(imm_does_not_pay_for_that_on_a_straight_target) {
    // The trade that matters. An IMM that only wins on turns has merely been
    // given more process noise; the point is to win there and not lose here.
    const Outcome o = run(0, 1000, 1000, 60, 991, false);
    std::printf("       straight run, no manoeuvre: IMM %.2f m, single-model EKF %.2f m\n",
                o.imm_rms, o.ekf_rms);
    std::printf("       peak P(manoeuvre) with nothing to detect: %.3f\n", o.peak_manoeuvre_prob);
    PT_CHECK(o.imm_rms < o.ekf_rms * 1.10);
    // A false manoeuvre alarm is cheap but should not be constant.
    PT_CHECK(o.peak_manoeuvre_prob < 0.9);
}

PT_TEST(imm_reports_which_way_the_target_turned) {
    // The sign of the turn-rate estimate is information a single-model filter
    // cannot produce at all.
    for (int dir = -1; dir <= 1; dir += 2) {
        TrackerConfig cfg;
        cfg.range_sigma_m = 5;
        cfg.bearing_sigma_rad = deg2rad(static_cast<Real>(0.5));
        ImmConfig imm;
        pt::Rng rng(70000 + static_cast<std::uint64_t>(dir + 1));
        Truth truth{2000, 8000, -30, 10};
        const Real omega = deg2rad(static_cast<Real>(3 * dir));

        ImmTrack it;
        imm_initiate(it, measure(truth, 0, rng, cfg), cfg, imm);
        Real w_est = 0;
        for (std::size_t k = 1; k <= 30; ++k) {
            advance(truth, (k >= 10) ? omega : static_cast<Real>(0), kDt);
            const Measurement z = measure(truth, static_cast<Real>(k) * kDt, rng, cfg);
            imm_predict(it, kDt, cfg, imm);
            imm_update(it, z, cfg, imm);
            if (k == 30) w_est = imm_turn_rate_estimate(it, imm);
        }
        std::printf("       truth %+d deg/s -> estimate %+.2f deg/s\n",
                    3 * dir, static_cast<double>(rad2deg(w_est)));
        PT_CHECK((static_cast<double>(w_est) > 0) == (dir > 0));
        PT_CHECK(std::fabs(static_cast<double>(rad2deg(w_est))) > 1.0);
    }
}

PT_TEST(model_probabilities_stay_a_probability_distribution) {
    // An invariant that is cheap to check and catastrophic to lose: the
    // probabilities must sum to one and none may reach zero, or a model can
    // never recover however well it would fit later.
    TrackerConfig cfg;
    ImmConfig imm;
    pt::Rng rng(31337);
    Truth truth{1500, 6000, -20, 5};
    ImmTrack it;
    imm_initiate(it, measure(truth, 0, rng, cfg), cfg, imm);

    double worst_sum_error = 0;
    Real smallest = 1;
    for (std::size_t k = 1; k <= 200; ++k) {
        // Deliberately violent: alternating hard turns, so every model in turn
        // becomes the badly-fitting one.
        const Real w = deg2rad(static_cast<Real>((k / 7) % 2 == 0 ? 6 : -6));
        advance(truth, w, kDt);
        const Measurement z = measure(truth, static_cast<Real>(k) * kDt, rng, cfg);
        imm_predict(it, kDt, cfg, imm);
        PT_CHECK(imm_update(it, z, cfg, imm));

        Real sum = 0;
        for (std::size_t j = 0; j < kImmModels; ++j) {
            sum += it.model_prob[j];
            smallest = std::min(smallest, it.model_prob[j]);
            PT_CHECK(it.model_prob[j] >= 0);
        }
        worst_sum_error = std::max(worst_sum_error, std::fabs(static_cast<double>(sum) - 1.0));
    }
    std::printf("       200 scans of alternating 6 deg/s turns:\n");
    std::printf("       worst |sum(mu) - 1| = %.2e, smallest model probability %.2e\n",
                worst_sum_error, static_cast<double>(smallest));
    PT_CHECK(worst_sum_error < pt::tol(1e-12, 1e-5));
    PT_CHECK(smallest > 0);
}

PT_TEST(the_zero_turn_rate_limit_is_exactly_constant_velocity) {
    // sin(w dt)/w is 0/0 at w = 0, and evaluating it directly loses all
    // precision long before w gets there. The implementation switches to a
    // series; this checks the switch is seamless, by walking the turn rate down
    // through the threshold and requiring the prediction to converge to the
    // straight-line answer rather than jump at the crossover.
    TrackerConfig cfg;
    Measurement z;
    z.range_m = 5000;
    z.bearing_rad = deg2rad(static_cast<Real>(20));
    z.time_s = 0;

    double previous = 1e30;
    for (int e = 2; e <= 9; ++e) {
        ImmConfig imm;
        imm.turn_rate_rps = static_cast<Real>(std::pow(10.0, -e));
        imm.switch_probability = static_cast<Real>(1e-6);   // keep the models apart
        ImmTrack it;
        imm_initiate(it, z, cfg, imm);
        for (std::size_t j = 0; j < kImmModels; ++j) {
            it.model_state[j].vx = 12;
            it.model_state[j].vy = -7;
        }
        imm_predict(it, 10, cfg, imm);

        // Model 1 turns at +omega; as omega -> 0 it must approach the straight
        // line the CV model (0) produces.
        const double d = std::sqrt(
            std::pow(static_cast<double>(it.model_state[1].x - it.model_state[0].x), 2) +
            std::pow(static_cast<double>(it.model_state[1].y - it.model_state[0].y), 2));
        std::printf("       omega = 1e-%d rad/s: |CT(+w) - CV| after 10 s = %.3e m\n", e, d);
        // Monotone convergence, but only while the difference is still
        // representable. The states are ~5000 m and Real may be float, whose
        // ~7 significant digits put a floor near 5e-4 m on any difference
        // between two such numbers. Below that the two models are the same
        // number and the comparison measures rounding, not the series.
        const double floor_m = pt::tol(1e-9, 1e-3);
        if (previous > floor_m) PT_CHECK(d < previous);
        previous = d;
    }
    PT_CHECK(previous < pt::tol(1e-5, 1e-3));
}
