// SPDX-License-Identifier: Apache-2.0
//
// Traces a Bellhop .env file with libphantom-sonar and writes the ray paths as
// CSV, so the same input file can be fed to both codes.
//
// Reading Bellhop's own input format rather than a transcribed copy is the
// point: it removes any possibility that the two tracers were shown different
// oceans. A transcription bug would look exactly like agreement.
//
//   phantom_trace <file.env> <out.csv> [--source-index N]
//
// This is tooling, not library code, so it uses std::string and std::vector
// freely. The zero-allocation guarantee covers include/ and src/.
#include "phantom/channel.hpp"
#include "phantom/profile.hpp"
#include "phantom/ray_tracer.hpp"

#include <array>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <utility>
#include <sstream>
#include <string>
#include <vector>

using namespace phantom;

namespace {

constexpr std::size_t kMaxProfilePoints = 16384;
constexpr std::size_t kMaxRayPoints     = 262144;

SoundSpeedProfile<kMaxProfilePoints> g_profile;
std::array<RayPoint, kMaxRayPoints>  g_scratch;

// --- minimal Fortran list-directed reader --------------------------------
//
// Handles the subset of the .env format the comparison needs: '!' comments,
// quoted strings, and '/' record terminators. Anything outside that subset is
// reported rather than guessed at.
class EnvReader {
  public:
    explicit EnvReader(const std::string& path) : in_(path) {}

    [[nodiscard]] bool ok() const { return in_.good() || !pending_.empty(); }

    // Next logical line with comments stripped.
    bool next_line(std::string& out) {
        std::string raw;
        while (std::getline(in_, raw)) {
            ++line_no_;
            const std::size_t bang = raw.find('!');
            if (bang != std::string::npos) raw.erase(bang);
            // Trim.
            const std::size_t b = raw.find_first_not_of(" \t\r\n");
            if (b == std::string::npos) continue;
            const std::size_t e = raw.find_last_not_of(" \t\r\n");
            out = raw.substr(b, e - b + 1);
            if (!out.empty()) return true;
        }
        return false;
    }

    [[nodiscard]] int line_no() const { return line_no_; }

  private:
    std::ifstream in_;
    std::string   pending_;
    int           line_no_ = 0;
};

// Strips a leading/trailing quote pair.
std::string unquote(std::string s) {
    if (s.size() >= 2 && (s.front() == '\'' || s.front() == '"')) {
        const char q = s.front();
        const std::size_t end = s.find(q, 1);
        if (end != std::string::npos) return s.substr(1, end - 1);
    }
    // Unquoted token: take up to the first space.
    const std::size_t sp = s.find_first_of(" \t");
    return sp == std::string::npos ? s : s.substr(0, sp);
}

// Reads the numbers on a line, stopping at a '/' terminator.
std::vector<double> numbers(const std::string& line) {
    std::vector<double> out;
    std::istringstream ss(line);
    std::string tok;
    while (ss >> tok) {
        if (tok == "/") break;
        if (!tok.empty() && tok.front() == '/') break;
        try {
            std::size_t used = 0;
            const double v = std::stod(tok, &used);
            if (used == 0) break;
            out.push_back(v);
        } catch (...) {
            break;
        }
    }
    return out;
}

struct EnvCase {
    std::string title;
    std::string ssp_option;
    double bottom_depth_m = 0;
    std::vector<std::pair<double, double>> ssp;  // (depth_m, speed_mps)
    std::vector<double> source_depths_m;
    std::vector<double> beam_angles_deg;
    double rbox_km = 0;
    double zbox_m  = 0;
    double step_m  = 0;
};

// Builds the profile, optionally splitting each tabulated layer into `refine`
// sublayers.
//
// This does NOT change the ocean. Within a tabulated layer c(z) is a straight
// line, so a point inserted on that line reproduces the same c(z) exactly, and
// `ray_is_invariant_under_layer_refinement` in the test suite pins that the
// traced answer does not move. What it changes is OUTPUT RESOLUTION: the tracer
// emits one point per layer crossing, so a coarse table yields a ray sampled
// every few kilometres. Comparing that against Bellhop by linear interpolation
// would measure the chord-versus-arc error of the comparison itself -- tens of
// metres on this profile -- rather than any disagreement between the codes.
bool build_profile(const EnvCase& env, std::size_t refine, std::string& err) {
    g_profile.clear();
    if (refine < 1) refine = 1;

    const std::size_t needed = (env.ssp.size() - 1) * refine + 1;
    if (needed > kMaxProfilePoints) {
        err = "refinement would need " + std::to_string(needed) + " profile points (max "
            + std::to_string(kMaxProfilePoints) + "); lower --refine";
        return false;
    }

    for (std::size_t i = 0; i + 1 < env.ssp.size(); ++i) {
        const double z0 = env.ssp[i].first,  c0 = env.ssp[i].second;
        const double z1 = env.ssp[i + 1].first, c1 = env.ssp[i + 1].second;
        for (std::size_t k = 0; k < refine; ++k) {
            const double f = static_cast<double>(k) / static_cast<double>(refine);
            if (!g_profile.push(static_cast<Real>(z0 + (z1 - z0) * f),
                                static_cast<Real>(c0 + (c1 - c0) * f))) {
                err = "SSP rejected at depth " + std::to_string(z0 + (z1 - z0) * f)
                    + " (non-monotonic depth or bad speed)";
                return false;
            }
        }
    }
    const auto& last = env.ssp.back();
    if (!g_profile.push(static_cast<Real>(last.first), static_cast<Real>(last.second))) {
        err = "SSP rejected at the bottom point";
        return false;
    }
    return true;
}

// Expands a Fortran "first last /" shorthand into `n` evenly spaced values.
std::vector<double> expand(const std::vector<double>& given, std::size_t n) {
    if (n == 0) return {};
    if (given.size() >= n) return std::vector<double>(given.begin(), given.begin() + static_cast<long>(n));
    if (given.size() == 1 || n == 1) return std::vector<double>(n, given.empty() ? 0.0 : given[0]);
    if (given.size() == 2) {
        std::vector<double> out(n);
        const double step = (given[1] - given[0]) / static_cast<double>(n - 1);
        for (std::size_t i = 0; i < n; ++i) out[i] = given[0] + step * static_cast<double>(i);
        out.back() = given[1];
        return out;
    }
    return given;
}

bool parse_env(const std::string& path, EnvCase& out, std::string& err) {
    EnvReader rd(path);
    if (!rd.ok()) { err = "cannot open " + path; return false; }

    std::string line;
    auto need = [&](const char* what) -> bool {
        if (!rd.next_line(line)) { err = std::string("unexpected end of file, wanted ") + what; return false; }
        return true;
    };

    if (!need("TITLE")) return false;
    out.title = unquote(line);
    if (!need("FREQ")) return false;
    if (!need("NMEDIA")) return false;
    if (!need("SSPOPT")) return false;
    out.ssp_option = unquote(line);

    if (!need("SSP header")) return false;
    const std::vector<double> hdr = numbers(line);
    if (hdr.size() < 3) { err = "malformed SSP header line"; return false; }
    out.bottom_depth_m = hdr[2];

    // SSP table, terminated when the tabulated depth reaches the bottom.
    out.ssp.clear();
    for (;;) {
        if (!need("SSP entry")) return false;
        const std::vector<double> v = numbers(line);
        if (v.size() < 2) { err = "malformed SSP entry: " + line; return false; }
        out.ssp.emplace_back(v[0], v[1]);
        if (v[0] >= out.bottom_depth_m - 1e-9) break;
    }
    if (out.ssp.size() < 2) { err = "SSP needs at least two points"; return false; }

    // Bottom boundary condition, plus its halfspace line when present.
    if (!need("bottom BC")) return false;
    const std::string bc = unquote(line);
    if (!bc.empty() && (bc[0] == 'A' || bc[0] == 'F' || bc[0] == 'G' || bc[0] == 'P')) {
        if (!need("halfspace parameters")) return false;
    }

    // Source depths.
    if (!need("NSD")) return false;
    const auto nsd = static_cast<std::size_t>(numbers(line).empty() ? 0 : numbers(line)[0]);
    if (!need("SD list")) return false;
    out.source_depths_m = expand(numbers(line), nsd);

    // Receiver depths and ranges are read only to advance the cursor: the
    // comparison is on ray geometry, not on a transmission-loss field.
    if (!need("NRD")) return false;
    if (!need("RD list")) return false;
    if (!need("NR")) return false;
    if (!need("R list")) return false;

    if (!need("run type")) return false;
    const std::string run_type = unquote(line);

    if (!need("NBEAMS")) return false;
    const std::vector<double> nb = numbers(line);
    const auto nbeams = static_cast<std::size_t>(nb.empty() ? 0 : nb[0]);
    if (nbeams == 0) {
        err = "NBEAMS = 0 means Bellhop picks the beam count itself, which "
              "would make the two runs incomparable. Set it explicitly.";
        return false;
    }

    if (!need("beam angles")) return false;
    out.beam_angles_deg = expand(numbers(line), nbeams);

    if (!need("STEP/ZBOX/RBOX")) return false;
    const std::vector<double> box = numbers(line);
    if (box.size() < 3) { err = "malformed STEP/ZBOX/RBOX line"; return false; }
    out.step_m  = box[0];
    out.zbox_m  = box[1];
    out.rbox_km = box[2];

    if (!run_type.empty() && run_type[0] != 'R' && run_type[0] != 'E') {
        std::fprintf(stderr,
                     "warning: run type '%s' is not a ray/eigenray run; "
                     "Bellhop will not write a .ray file\n",
                     run_type.c_str());
    }
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 3) {
        std::fprintf(stderr,
                     "usage: %s <file.env> <out.csv> [--source-index N] [--refine N]\n"
                     "\n"
                     "  --refine N   split each tabulated SSP layer into N sublayers.\n"
                     "               Does not change the modelled ocean (c(z) is linear\n"
                     "               within a layer); raises the ray sampling density so\n"
                     "               the polyline can be compared against a stepping code\n"
                     "               without chord-versus-arc error. Default 32.\n",
                     argv[0]);
        return 2;
    }
    const std::string env_path = argv[1];
    const std::string csv_path = argv[2];
    std::size_t source_index = 0;
    std::size_t refine = 32;
    for (int i = 3; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--source-index" && i + 1 < argc) {
            source_index = static_cast<std::size_t>(std::atol(argv[++i]));
        } else if (a == "--refine" && i + 1 < argc) {
            refine = static_cast<std::size_t>(std::atol(argv[++i]));
        }
    }

    EnvCase env;
    std::string err;
    if (!parse_env(env_path, env, err)) {
        std::fprintf(stderr, "error: %s\n", err.c_str());
        return 1;
    }

    // The interpolation scheme has to match or the two codes are modelling
    // different oceans. libphantom-sonar is piecewise linear in c(z), which is
    // Bellhop's 'C' option. PCHIP ('P') and spline ('S') are smooth
    // reconstructions of the same table and will legitimately disagree.
    const char interp = env.ssp_option.empty() ? '?' : env.ssp_option[0];
    if (interp != 'C') {
        std::fprintf(stderr,
                     "error: SSP interpolation is '%c'; libphantom-sonar is "
                     "piecewise linear, which is Bellhop's 'C'. Comparing "
                     "against '%c' would measure the interpolation difference, "
                     "not the tracer.\n", interp, interp);
        return 1;
    }

    if (env.source_depths_m.empty()) {
        std::fprintf(stderr, "error: no source depths in %s\n", env_path.c_str());
        return 1;
    }
    if (source_index >= env.source_depths_m.size()) {
        std::fprintf(stderr, "error: --source-index %zu but only %zu sources\n",
                     source_index, env.source_depths_m.size());
        return 1;
    }

    if (!build_profile(env, refine, err)) {
        std::fprintf(stderr, "error: %s\n", err.c_str());
        return 1;
    }

    const ProfileView svp = g_profile.view();
    const auto source_depth = static_cast<Real>(env.source_depths_m[source_index]);

    TraceConfig cfg;
    cfg.max_range_m    = static_cast<Real>(env.rbox_km * 1000.0);
    cfg.max_time_s     = static_cast<Real>(1e7);   // range is the binding limit
    cfg.bottom_depth_m = static_cast<Real>(env.bottom_depth_m);
    cfg.surface        = BoundaryAction::Reflect;
    cfg.bottom         = BoundaryAction::Reflect;
    cfg.max_bounces    = 100000;

    std::FILE* f = std::fopen(csv_path.c_str(), "w");
    if (f == nullptr) {
        std::fprintf(stderr, "error: cannot write %s\n", csv_path.c_str());
        return 1;
    }
    std::fprintf(f, "ray_id,launch_deg,range_m,depth_m,angle_deg,time_s,speed_mps\n");

    std::size_t total_points = 0;
    std::size_t total_steps  = 0;
    std::size_t truncated    = 0;

    for (std::size_t i = 0; i < env.beam_angles_deg.size(); ++i) {
        const auto launch = static_cast<Real>(env.beam_angles_deg[i]);
        const TraceResult r =
            trace_ray(svp, source_depth, deg2rad(launch), cfg, g_scratch);
        if (r.status == TraceStatus::BufferFull) ++truncated;
        total_points += r.point_count;
        total_steps  += r.steps;
        for (std::size_t k = 0; k < r.point_count; ++k) {
            const RayPoint& p = g_scratch[k];
            std::fprintf(f, "%zu,%.6f,%.6f,%.6f,%.6f,%.9f,%.6f\n",
                         i, static_cast<double>(launch),
                         static_cast<double>(p.range_m),
                         static_cast<double>(p.depth_m),
                         static_cast<double>(rad2deg(p.angle_rad)),
                         static_cast<double>(p.time_s),
                         static_cast<double>(p.speed_mps));
        }
    }
    std::fclose(f);

    const ChannelInfo ch = analyze_sofar(svp);

    std::printf("phantom_trace: %s\n", env.title.c_str());
    std::printf("  env file       : %s\n", env_path.c_str());
    std::printf("  SSP            : %zu tabulated points, 0 - %.1f m, interpolation '%c'\n",
                env.ssp.size(), static_cast<double>(env.bottom_depth_m), interp);
    std::printf("  refinement     : x%zu -> %zu profile points (same c(z), denser output)\n",
                refine, svp.point_count());
    if (ch.found) {
        std::printf("  channel axis   : %.1f m at %.2f m/s, cone +/- %.3f deg\n",
                    static_cast<double>(ch.axis_depth_m),
                    static_cast<double>(ch.axis_speed_mps),
                    static_cast<double>(rad2deg(ch.max_trapped_angle_rad)));
    }
    std::printf("  source depth   : %.1f m (index %zu of %zu)\n",
                env.source_depths_m[source_index], source_index,
                env.source_depths_m.size());
    std::printf("  beams          : %zu from %.3f to %.3f deg\n",
                env.beam_angles_deg.size(), env.beam_angles_deg.front(),
                env.beam_angles_deg.back());
    std::printf("  range limit    : %.1f km\n", env.rbox_km);
    std::printf("  arc steps      : %zu, polyline points %zu\n", total_steps, total_points);
    if (truncated > 0) {
        std::fprintf(stderr, "warning: %zu ray(s) hit the point-buffer limit\n", truncated);
    }
    std::printf("  wrote          : %s\n", csv_path.c_str());
    return 0;
}
