// RiskIndex: truncated-spread risk basis  F_tau(S) = (1/M) sum_m min{ hat X_m(S), tau }
// and threshold-specific lazy-greedy chains C_h = (u_{h,1}, ..., u_{h,Kmax}).
#pragma once
#include "../grr/grr.hpp"
#include <queue>
#include <limits>

namespace ri {

static const double TAU_INF = std::numeric_limits<double>::infinity();

// Coverage state of one greedy chain over a GRR sketch.
struct CoverState {
    const GRRSketch& sk;
    uint32_t MB;                    // number of *build* worlds (groups 0..MB-1) used for gains
    std::vector<uint32_t> cnt;      // c_m = #covered RR sets in group m
    std::vector<uint8_t> covered;   // per RR id
    std::vector<uint32_t> delta;    // scratch: delta_m(u | S)
    std::vector<uint32_t> touched;  // scratch: groups touched by u
    explicit CoverState(const GRRSketch& s, uint32_t M_build = 0) : sk(s), MB(M_build ? M_build : s.M), cnt(s.M, 0), covered(s.num_sets(), 0), delta(s.M, 0) {}

    // delta_m(u|S) for all m into `delta`/`touched`; returns total uncovered incidences
    inline uint32_t compute_delta(uint32_t u) {
        touched.clear(); uint32_t tot = 0;
        for (uint64_t i = sk.inv_off[u]; i < sk.inv_off[u + 1]; ++i) {
            uint32_t r = sk.inv_rr[i];
            if (covered[r]) continue;
            uint32_t m = r / sk.b;
            if (delta[m] == 0) touched.push_back(m);
            delta[m]++; ++tot;
        }
        return tot;
    }
    inline void clear_delta() { for (uint32_t m : touched) delta[m] = 0; touched.clear(); }

    // Eq. (marginal basis gain): (1/M) sum_m [ min{(n/b)(c_m+delta_m), tau} - min{(n/b)c_m, tau} ]
    // (restricted to the MB build worlds when cross-fitting is enabled)
    inline double gain_tau(double tau) const {
        double sc = sk.scale(), g = 0;
        if (tau == TAU_INF) { for (uint32_t m : touched) if (m < MB) g += (double)delta[m]; return g * sc / MB; }
        for (uint32_t m : touched) {
            if (m >= MB) continue;
            double before = std::min(sc * cnt[m], tau);
            double after = std::min(sc * (cnt[m] + delta[m]), tau);
            g += after - before;
        }
        return g / MB;
    }
    inline void commit(uint32_t u) {  // add u to S
        for (uint64_t i = sk.inv_off[u]; i < sk.inv_off[u + 1]; ++i) {
            uint32_t r = sk.inv_rr[i];
            if (!covered[r]) { covered[r] = 1; cnt[r / sk.b]++; }
        }
    }
    // profile over the *validation* worlds MB..M-1 (all worlds when MB == M)
    inline void profile(std::vector<double>& X) const {
        double sc = sk.scale(); uint32_t m0 = (MB == sk.M) ? 0 : MB; X.resize(sk.M - m0);
        for (uint32_t m = m0; m < sk.M; ++m) X[m - m0] = sc * cnt[m];
    }
};

struct ChainResult {
    double tau;
    std::vector<uint32_t> seeds;                 // u_{h,1..K}
    std::vector<std::vector<double>> prefix_sum; // [k-1] -> P_0..P_M of the sorted profile of S_h^{(k)}
    std::vector<std::vector<double>> prefix_sq;  // [k-1] -> Q_0..Q_M prefix sums of squares (optional; tail standard errors)
    std::vector<double> basis_value;             // hat F_tau(S_h^{(k)})
    std::vector<double> rho;                     // tail-resolution ratio of each prefix profile
    std::vector<double> mean;                    // sample mean of each prefix profile (preserved by variance matching)
    uint64_t gain_evals = 0;
};

// Materialize the sorted (optionally variance-matched) scenario profile of the
// current prefix held by `st` into `res`: prefix sums, prefix sums of squares,
// resolution ratio rho and profile mean.
inline void append_profile(ChainResult& res, const CoverState& st, bool vm, bool store_sq, std::vector<double>& X) {
    st.profile(X);
    res.mean.push_back(mean_of(X));
    res.rho.push_back(vm ? vm_shrink(X, (double)st.sk.n, (double)st.sk.b) : 1.0);
    std::sort(X.begin(), X.end());
    std::vector<double> P(X.size() + 1, 0.0);
    for (size_t m = 0; m < X.size(); ++m) P[m + 1] = P[m] + X[m];
    res.prefix_sum.push_back(std::move(P));
    if (store_sq) {
        std::vector<double> Q(X.size() + 1, 0.0);
        for (size_t m = 0; m < X.size(); ++m) Q[m + 1] = Q[m] + X[m] * X[m];
        res.prefix_sq.push_back(std::move(Q));
    }
}

// Lazy (CELF) greedy for the submodular sketched objective hat F_tau.
inline ChainResult build_chain(const GRRSketch& sk, double tau, uint32_t K, bool materialize = true, uint32_t M_build = 0, bool vm = true, bool store_sq = false) {
    ChainResult res; res.tau = tau;
    CoverState st(sk, M_build);
    struct Item { double ub; uint32_t u; uint32_t it; bool operator<(const Item& o) const { return ub < o.ub; } };
    std::priority_queue<Item> pq;
    // initial exact gains at S = empty
    for (uint32_t u = 0; u < sk.n; ++u) {
        if (sk.deg(u) == 0) continue;
        st.compute_delta(u); double g = st.gain_tau(tau); st.clear_delta(); res.gain_evals++;
        pq.push({g, u, 1});  // exact at S=empty, hence fresh for iteration 1
    }
    std::vector<uint8_t> chosen(sk.n, 0);
    std::vector<double> X;
    double value = 0;
    for (uint32_t k = 1; k <= K && !pq.empty(); ++k) {
        while (true) {
            Item top = pq.top(); pq.pop();
            if (top.it == k) {  // fresh evaluation -> select
                chosen[top.u] = 1; st.commit(top.u); res.seeds.push_back(top.u); value += top.ub;
                break;
            }
            st.compute_delta(top.u); double g = st.gain_tau(tau); st.clear_delta(); res.gain_evals++;
            pq.push({g, top.u, k});
        }
        res.basis_value.push_back(value);
        if (materialize) append_profile(res, st, vm, store_sq, X);
    }
    return res;
}

// Hybrid mean chain: a seed sequence chosen elsewhere (standard RIS greedy on an
// independent RR sketch, which pays no price of grouping for the expected spread)
// is materialized on the GRR sketch, so that its prefix profiles and certificates
// are comparable with those of the threshold chains.
inline ChainResult materialize_sequence(const GRRSketch& sk, const std::vector<uint32_t>& seq, bool vm = true, bool store_sq = false) {
    ChainResult res; res.tau = TAU_INF;
    CoverState st(sk);
    std::vector<double> X;
    double value = 0;
    for (uint32_t u : seq) {
        st.compute_delta(u); value += st.gain_tau(TAU_INF); st.clear_delta();
        st.commit(u); res.seeds.push_back(u); res.basis_value.push_back(value);
        append_profile(res, st, vm, store_sq, X);
    }
    return res;
}

// Standard error of the certificate including the within-world noise of the
// tail scenarios: se_w^2 = se^2 + rho^2 * sum_{tail} v_m / r^2, where
// v_m = x_m (n - x_m) / (b - 1) is the unbiased within-world variance of the raw
// scenario estimate x_m (Prop. vm) and x_m is recovered from the stored corrected
// value z_m = mu + rho (x_m - mu).  O(1) from (P, Q, mu, rho, n, b).
static inline double tail_se_within(const std::vector<double>& P, const std::vector<double>& Q, double alpha,
                                    double mu, double rho, double n, double b) {
    size_t M = P.size() - 1;
    double aM = alpha * (double)M;
    size_t r = (size_t)std::ceil(aM - 1e-12); if (r == 0) r = 1; if (r > M) r = M;
    double sz = P[r], szz = Q[r], rr = (double)r;
    double mean = sz / rr, var = szz / rr - mean * mean; if (var < 0) var = 0;
    double se2 = var / rr;
    if (rho > 0 && b > 1) {
        double sx = rr * mu + (sz - rr * mu) / rho;
        double sxx = rr * mu * mu + 2.0 * mu * (sz - rr * mu) / rho + (szz - 2.0 * mu * sz + rr * mu * mu) / (rho * rho);
        double sv = (n * sx - sxx) / (b - 1.0); if (sv < 0) sv = 0;
        se2 += rho * rho * sv / (rr * rr);
    }
    return std::sqrt(se2);
}

// Standard error of the empirical lower-tail mean (the r = ceil(alpha M) lowest
// scenarios) from prefix sums P and prefix sums of squares Q, in O(1):
// sqrt( Var(z_(1..r)) / r ).  A plug-in dispersion measure of the certificate's
// between-world sampling noise (it ignores the VaR-boundary term and treats the
// tail members as fixed), used only as a scale for the selection margin and as a
// reported uncertainty; not a confidence bound.
static inline double tail_se_from_prefix(const std::vector<double>& P, const std::vector<double>& Q, double alpha) {
    size_t M = P.size() - 1;
    double aM = alpha * (double)M;
    size_t r = (size_t)std::ceil(aM - 1e-12); if (r == 0) r = 1; if (r > M) r = M;
    double mean = P[r] / r, var = Q[r] / r - mean * mean; if (var < 0) var = 0;
    return std::sqrt(var / (double)r);
}

// Empirical lower-tail CVaR from a prefix-sum profile in O(1).
static inline double cvar_from_prefix(const std::vector<double>& P, double alpha) {
    size_t M = P.size() - 1;
    double aM = alpha * (double)M;
    size_t r = (size_t)std::ceil(aM - 1e-12); if (r == 0) r = 1; if (r > M) r = M;
    double z_r = P[r] - P[r - 1];
    return (P[r - 1] + (aM - (double)(r - 1)) * z_r) / aM;
}

}  // namespace ri
