// RiskIndex: common utilities (RNG, hashing, timing)
#pragma once
#include <cstdint>
#include <chrono>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <sys/resource.h>

namespace ri {

// --- deterministic 64-bit mixing (splitmix64 finalizer) ----------------------
static inline uint64_t mix64(uint64_t z) {
    z += 0x9E3779B97F4A7C15ULL;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}
// uniform double in [0,1) from (stream, index) pair -- used for scenario-consistent
// edge/node states: the state of edge e in world w is a pure function of (w, e).
static inline double hash_uniform(uint64_t stream, uint64_t idx) {
    uint64_t h = mix64(stream ^ mix64(idx + 0x1234567ULL));
    return (h >> 11) * (1.0 / 9007199254740992.0);
}

// --- xoshiro-like PCG RNG for sequential sampling ----------------------------
struct Rng {
    uint64_t s;
    explicit Rng(uint64_t seed) : s(mix64(seed ^ 0xA5A5A5A5A5A5A5A5ULL)) {}
    inline uint64_t next() { s = mix64(s); return s; }
    inline double uniform() { return (next() >> 11) * (1.0 / 9007199254740992.0); }
    inline uint32_t below(uint32_t n) { return (uint32_t)(uniform() * n); }
};

// --- timing ------------------------------------------------------------------
struct Timer {
    std::chrono::steady_clock::time_point t0;
    Timer() : t0(std::chrono::steady_clock::now()) {}
    double sec() const {
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    }
    void reset() { t0 = std::chrono::steady_clock::now(); }
};

static inline double peak_rss_mb() {
    struct rusage ru; getrusage(RUSAGE_SELF, &ru);
    return ru.ru_maxrss / 1024.0;  // Linux: KB -> MB
}

// --- lower-tail CVaR of an empirical sample ----------------------------------
// Exact empirical lower-tail CVaR at level alpha (Rockafellar-Uryasev convention
// for rewards): average of the lowest alpha-fraction of the mass, with the
// fractional correction at the boundary order statistic.  Requires sorted z.
static inline double cvar_sorted(const std::vector<double>& z, double alpha) {
    size_t M = z.size();
    double aM = alpha * (double)M;
    size_t r = (size_t)std::ceil(aM - 1e-12);
    if (r == 0) r = 1;
    if (r > M) r = M;
    double s = 0;
    for (size_t i = 0; i + 1 < r; ++i) s += z[i];
    s += (aM - (double)(r - 1)) * z[r - 1];
    return s / aM;
}
static inline double cvar_unsorted(std::vector<double> z, double alpha) {
    std::sort(z.begin(), z.end());
    return cvar_sorted(z, alpha);
}
// --- variance-matched (noise-aware) scenario profile ---------------------------
// Each hat X_m = X_m + e_m with known within-world sampling variance
// Var(e_m | L_m) = X_m (n - X_m) / b  (b targets with replacement).  The unbiased
// plug-in estimate v_m = hat X_m (n - hat X_m)/(b-1) satisfies E[v_m | L_m] =
// X_m (n - X_m)/b, so  rho^2 = 1 - mean(v_m)/S^2  estimates Var(X)/Var(hat X).
// The shrunk profile  mu + rho (hat X_m - mu)  keeps the mean and matches the
// variance of the true world distribution, removing the downward tail bias that
// independent per-world noise induces in the lowest order statistics.
// Returns rho in [0,1] (0 = no resolvable scenario variation at this b).
static inline double vm_shrink(std::vector<double>& X, double n, double b) {
    size_t M = X.size(); if (M < 2 || b <= 1) return 1.0;
    double mu = 0; for (double x : X) mu += x; mu /= M;
    double S2 = 0, vbar = 0;
    for (double x : X) { S2 += (x - mu) * (x - mu); double xc = std::min(std::max(x, 0.0), n); vbar += xc * (n - xc) / (b - 1.0); }
    S2 /= (M - 1); vbar /= M;
    double rho = (S2 > vbar && S2 > 0) ? std::sqrt(1.0 - vbar / S2) : 0.0;
    for (double& x : X) x = mu + rho * (x - mu);
    return rho;
}
static inline double mean_of(const std::vector<double>& z) {
    double s = 0; for (double v : z) s += v; return z.empty() ? 0 : s / z.size();
}

}  // namespace ri
