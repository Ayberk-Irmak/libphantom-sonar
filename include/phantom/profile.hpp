// SPDX-License-Identifier: Apache-2.0
// libphantom-sonar — stratified sound speed profile with precomputed gradients.
//
// The profile is stored as a piecewise-LINEAR function of depth, not as a stack
// of isovelocity slabs. That single choice is what lets the ray tracer use exact
// circular arcs instead of stepping, which is both faster and more accurate.
#ifndef PHANTOM_PROFILE_HPP
#define PHANTOM_PROFILE_HPP

#include "phantom/types.hpp"

#include <array>
#include <span>

namespace phantom {

// Non-owning view over a finalized profile. This is what the hot path consumes,
// so the kernel compiles once regardless of the caller's capacity template.
struct ProfileView {
    std::span<const Real> depth_m;    // strictly increasing, size n >= 2
    std::span<const Real> speed_mps;  // size n
    std::span<const Real> gradient;   // size n - 1, dc/dz within each layer

    [[nodiscard]] constexpr std::size_t point_count() const noexcept { return depth_m.size(); }
    [[nodiscard]] constexpr std::size_t layer_count() const noexcept { return gradient.size(); }
    [[nodiscard]] constexpr Real min_depth() const noexcept { return depth_m.front(); }
    [[nodiscard]] constexpr Real max_depth() const noexcept { return depth_m.back(); }
    [[nodiscard]] constexpr bool valid() const noexcept {
        return depth_m.size() >= 2 && speed_mps.size() == depth_m.size()
            && gradient.size() + 1 == depth_m.size();
    }
};

// Index of the layer [depth[i], depth[i+1]] that contains `z`. When `z` sits
// exactly on a node, `downgoing` selects which side to take -- getting this
// wrong is how ray tracers stall on layer boundaries.
[[nodiscard]] inline std::size_t find_layer(const ProfileView& svp, Real z, bool downgoing) noexcept {
    const std::size_t n = svp.depth_m.size();
    std::size_t lo = 0;
    std::size_t hi = n - 1;
    while (lo + 1 < hi) {
        const std::size_t mid = lo + (hi - lo) / 2;
        if (svp.depth_m[mid] <= z) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    if (!downgoing && lo > 0 && svp.depth_m[lo] >= z) {
        --lo;
    }
    const std::size_t last = n - 2;
    return lo > last ? last : lo;
}

// Linear interpolation of sound speed. Clamps outside the profile.
[[nodiscard]] inline Real speed_at(const ProfileView& svp, Real z) noexcept {
    const Real zc = clamp(z, svp.min_depth(), svp.max_depth());
    const std::size_t i = find_layer(svp, zc, true);
    return svp.speed_mps[i] + svp.gradient[i] * (zc - svp.depth_m[i]);
}

// Owning, fixed-capacity profile. No heap, ever.
template <std::size_t Capacity>
class SoundSpeedProfile {
    static_assert(Capacity >= 2, "a profile needs at least two points");

  public:
    // Appends a sample. Depths must be strictly increasing. Returns false on
    // capacity overflow or on a non-monotonic / non-physical sample.
    bool push(Real depth_m, Real speed_mps) noexcept {
        if (count_ >= Capacity) return false;
        if (!(speed_mps > 0)) return false;
        if (count_ > 0 && !(depth_m > depth_[count_ - 1])) return false;
        depth_[count_] = depth_m;
        speed_[count_] = speed_mps;
        if (count_ > 0) {
            const std::size_t seg = count_ - 1;
            grad_[seg] = (speed_[count_] - speed_[seg]) / (depth_[count_] - depth_[seg]);
        }
        ++count_;
        return true;
    }

    void clear() noexcept { count_ = 0; }

    [[nodiscard]] std::size_t size() const noexcept { return count_; }
    [[nodiscard]] bool valid() const noexcept { return count_ >= 2; }
    [[nodiscard]] static constexpr std::size_t capacity() noexcept { return Capacity; }

    [[nodiscard]] ProfileView view() const noexcept {
        return ProfileView{
            std::span<const Real>(depth_.data(), count_),
            std::span<const Real>(speed_.data(), count_),
            std::span<const Real>(grad_.data(), count_ > 0 ? count_ - 1 : 0),
        };
    }

  private:
    std::array<Real, Capacity>     depth_{};
    std::array<Real, Capacity>     speed_{};
    std::array<Real, Capacity - 1> grad_{};
    std::size_t                    count_ = 0;
};

// Samples an analytic c(z) onto a profile at `n` evenly spaced depths.
// Useful for the Munk benchmark and for synthetic test cases.
template <std::size_t Capacity, typename SpeedFn>
bool fill_profile(SoundSpeedProfile<Capacity>& out, Real z_top, Real z_bottom,
                  std::size_t n, SpeedFn&& c_of_z) {
    if (n < 2 || n > Capacity || !(z_bottom > z_top)) return false;
    out.clear();
    const Real dz = (z_bottom - z_top) / static_cast<Real>(n - 1);
    for (std::size_t i = 0; i < n; ++i) {
        const Real z = (i + 1 == n) ? z_bottom : z_top + dz * static_cast<Real>(i);
        if (!out.push(z, c_of_z(z))) return false;
    }
    return true;
}

}  // namespace phantom

#endif  // PHANTOM_PROFILE_HPP
