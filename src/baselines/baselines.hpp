// RiskIndex baselines operating on a (grouped or independent) RR sketch.
//
//  * direct_cvar_greedy : per-query greedy on the empirical lower-tail CVaR of the
//    scenario profile (non-submodular).  Lazy evaluation uses a *valid* upper
//    bound on the CVaR gain, with stale (submodular, hence admissible) mean gains:
//      raw objective:  DeltaCVaR(u|S) <= (1/alpha) DeltaMean(u|S)
//        (CVaR_alpha(Y) = min over tail weightings w, w_m <= 1/(alpha M), of
//         sum_m w_m Y_m; evaluating the optimal weighting of X on Y = X + d gives
//         CVaR(Y) <= CVaR(X) + (1/alpha) mean(d)).
//      variance-matched objective:  the corrected CVaR of any profile Y is
//        mu_Y + rho_Y (CVaR_raw(Y) - mu_Y) <= mu_Y, because rho_Y in [0,1] and
//        CVaR_raw <= mu.  Hence DeltaCVaR_VM(u|S) <= [mu_S - CVaR_VM(S)] + DeltaMean(u|S).
//        The (1/alpha)-bound is NOT valid here: rho_Y depends on the candidate,
//        and a candidate that raises rho can raise the corrected CVaR by more
//        than (1/alpha) times its mean gain.  (Earlier versions used the raw bound
//        for the VM objective; see experiments/check_vm_pruning for the audit.)
//    At alpha = 1 the raw variant reduces exactly to lazy greedy on expected
//    coverage.  The CVaR of a candidate profile is computed by merging the changed
//    groups into the sorted current profile, O(r + t log t).
//  * mean_greedy       : the classical RIS / IMM-style selection (tau = inf chain).
//  * seq_greedy        : the sequential greedy algorithm (SGA) of Zhou & Tokekar
//    [zhou2020] on the same sketch: a uniform threshold grid tau_i = i*Gamma/H',
//    one greedy run per threshold on the surrogate H(S,tau) = tau - (1/alpha) E[(tau - X)^+]
//    (equivalently (1/alpha) F_tau(S) - (1-alpha)/alpha tau, which is submodular
//    in S for fixed tau), and the pair (S, tau) with the largest surrogate value
//    is returned.  The surrogate is evaluated on the M worlds of the sketch, which
//    play the role of the Monte Carlo samples of [zhou2020].  Gamma is the largest
//    per-world spread of the mean-greedy k-set on the sketch (the surrogate is
//    constant for tau beyond the support).  select_cvar=true selects by the exact
//    empirical CVaR of the candidate instead (a strengthened variant).
#pragma once
#include "../risk_basis/risk_basis.hpp"

namespace ri {

struct GreedyStats { uint64_t cvar_evals = 0, mean_evals = 0; double sec = 0; };

inline std::vector<uint32_t> direct_cvar_greedy(const GRRSketch& sk, uint32_t k, double alpha, GreedyStats* stats = nullptr, bool vm = false,
                                                bool legacy_vm_bound = false, bool no_prune = false) {
    Timer T;
    CoverState st(sk);
    const double sc = sk.scale(), inv_alpha = 1.0 / alpha;
    const uint32_t M = sk.M;
    const double aM = alpha * (double)M;
    size_t r = (size_t)std::ceil(aM - 1e-12); if (r == 0) r = 1; if (r > M) r = M;
    // Item.ub holds the (stale) mean gain; the admissible CVaR-gain bound is
    // derived from it on comparison (bound(mg) below).
    struct Item { double ub; uint32_t u; uint32_t it; bool operator<(const Item& o) const { return ub < o.ub; } };
    std::priority_queue<Item> pq;
    uint64_t cvar_evals = 0, mean_evals = 0;
    for (uint32_t u = 0; u < sk.n; ++u) if (sk.deg(u) > 0) pq.push({sc * sk.deg(u) / M, u, 0});
    std::vector<double> cur(M, 0.0);
    std::vector<std::pair<double, uint32_t>> sorted(M);  // (value, group) ascending
    for (uint32_t m = 0; m < M; ++m) sorted[m] = {0.0, m};
    std::vector<uint8_t> is_touched(M, 0);
    std::vector<double> newvals;
    double cur_cvar = 0.0, cur_mean = 0.0;
    std::vector<uint32_t> S;
    std::vector<double> full;
    auto tail_cvar = [&](const std::vector<uint32_t>& touched) {
        if (vm) {  // variance-matched objective: full recompute of the candidate profile
            full = cur; for (uint32_t m : touched) full[m] += sc * st.delta[m];
            vm_shrink(full, (double)sk.n, (double)sk.b); return cvar_unsorted(full, alpha);
        }
        newvals.clear();
        for (uint32_t m : touched) { is_touched[m] = 1; newvals.push_back(cur[m] + sc * st.delta[m]); }
        std::sort(newvals.begin(), newvals.end());
        size_t i = 0, j = 0, taken = 0; double s = 0, last = 0;
        while (taken < r) {
            while (i < M && is_touched[sorted[i].second]) ++i;
            double v;
            if (j < newvals.size() && (i >= M || newvals[j] <= sorted[i].first)) v = newvals[j++];
            else v = sorted[i++].first;
            ++taken;
            if (taken < r) s += v; else last = v;
        }
        for (uint32_t m : touched) is_touched[m] = 0;
        return (s + (aM - (double)(r - 1)) * last) / aM;
    };
    for (uint32_t k_ = 1; k_ <= k && !pq.empty(); ++k_) {
        const double gap = cur_mean - cur_cvar;  // >= 0
        auto bound = [&](double mg) { if (no_prune) return std::numeric_limits<double>::infinity(); return (vm && !legacy_vm_bound) ? mg + gap : inv_alpha * mg; };
        double best_gain = -1; uint32_t best_u = UINT32_MAX;
        std::vector<Item> stale;
        while (!pq.empty()) {
            Item top = pq.top();
            if (bound(top.ub) <= best_gain) break;  // all remaining true CVaR gains <= best_gain
            pq.pop();
            uint32_t tot = st.compute_delta(top.u); ++mean_evals;
            double mean_gain = sc * (double)tot / M;
            double c = tail_cvar(st.touched); ++cvar_evals;
            st.clear_delta();
            double g = c - cur_cvar;
            if (g > best_gain) { best_gain = g; best_u = top.u; }
            stale.push_back({mean_gain, top.u, k_});
        }
        if (best_u == UINT32_MAX) break;
        for (auto& it : stale) if (it.u != best_u) pq.push(it);
        st.compute_delta(best_u);
        for (uint32_t m : st.touched) cur[m] += sc * st.delta[m];
        st.clear_delta();
        st.commit(best_u);
        for (uint32_t m = 0; m < M; ++m) sorted[m] = {cur[m], m};
        std::sort(sorted.begin(), sorted.end());
        cur_cvar = 0; { double s = 0; for (size_t i = 0; i + 1 < r; ++i) s += sorted[i].first; cur_cvar = (s + (aM - (double)(r - 1)) * sorted[r - 1].first) / aM; }
        if (vm) { full = cur; vm_shrink(full, (double)sk.n, (double)sk.b); cur_cvar = cvar_unsorted(full, alpha); }
        cur_mean = mean_of(cur);
        S.push_back(best_u);
    }
    if (stats) { stats->cvar_evals = cvar_evals; stats->mean_evals = mean_evals; stats->sec = T.sec(); }
    return S;
}

inline std::vector<uint32_t> mean_greedy(const GRRSketch& sk, uint32_t k) {
    ChainResult c = build_chain(sk, TAU_INF, k, false);
    return c.seeds;
}

// Surrogate value H(S,tau) = (1/alpha) F_tau(S) - (1-alpha)/alpha * tau of a profile X.
static inline double surrogate_value(const std::vector<double>& X, double tau, double alpha) {
    if (tau == TAU_INF) return mean_of(X);  // alpha = 1 member; for alpha < 1 the surrogate is -inf at tau = inf
    double s = 0; for (double x : X) s += std::min(x, tau);
    return (s / X.size()) / alpha - (1.0 - alpha) / alpha * tau;
}

// Uniform grid of [zhou2020]: tau_i = i * Gamma / Hgrid, i = 1..Hgrid.
static inline std::vector<double> uniform_grid(double Gamma, uint32_t Hgrid) {
    std::vector<double> taus;
    for (uint32_t i = 1; i <= Hgrid; ++i) taus.push_back(Gamma * (double)i / (double)Hgrid);
    return taus;
}

struct SeqGreedyResult { std::vector<uint32_t> S; double tau = 0; double value = 0; GreedyStats stats; };

// Sequential greedy [zhou2020] for one query (k, alpha): Hgrid greedy runs to k.
inline SeqGreedyResult seq_greedy(const GRRSketch& sk, uint32_t k, double alpha, uint32_t Hgrid, bool select_cvar = false,
                                  const std::vector<uint32_t>* mean_set = nullptr) {
    Timer T; SeqGreedyResult res;
    std::vector<double> X; std::vector<uint8_t> covered(sk.num_sets(), 0); std::vector<uint32_t> touched;
    std::vector<uint32_t> ms;
    if (mean_set) ms = *mean_set; else ms = mean_greedy(sk, k);
    sk.profile(ms, X, covered, touched);
    double Gamma = *std::max_element(X.begin(), X.end());
    if (Gamma <= 0) Gamma = sk.scale();
    auto taus = uniform_grid(Gamma, Hgrid);
    double best = -std::numeric_limits<double>::infinity();
    for (double tau : taus) {
        ChainResult c = build_chain(sk, tau, k, false);
        res.stats.mean_evals += c.gain_evals;
        sk.profile(c.seeds, X, covered, touched);
        double v = select_cvar ? cvar_unsorted(X, alpha) : surrogate_value(X, tau, alpha);
        if (v > best) { best = v; res.S = c.seeds; res.tau = tau; res.value = v; }
    }
    // alpha = 1 has no finite maximizer of the surrogate; [zhou2020] assumes alpha < 1.
    // We include the mean set as a candidate so that the baseline is well defined at alpha = 1.
    if (alpha >= 1.0) { res.S = ms; res.tau = TAU_INF; res.value = mean_of(X); }
    res.stats.sec = T.sec();
    return res;
}

// Cached sequential greedy: the Hgrid chains of the uniform grid are built once to
// Kmax (grid from the Kmax mean prefix, independent of the query) and every query
// (k, alpha) re-evaluates the surrogate of each k-prefix from the sketch.  This is
// the natural "reuse" version of [zhou2020]; it differs from RiskIndex in grid
// (uniform vs. sketch-derived geometric), in selection (surrogate vs. exact CVaR
// of the variance-matched profile) and in materialization (none).
struct SeqGreedyCache {
    const GRRSketch& sk;
    std::vector<double> taus;
    std::vector<ChainResult> chains;
    std::vector<uint32_t> mean_chain;
    double build_sec = 0;
    std::vector<uint8_t> covered; std::vector<uint32_t> touched; std::vector<double> X;
    // grid == nullptr: uniform grid of [zhou2020] on (0, Gamma]; otherwise the given thresholds (finite ones;
    // tau = inf, i.e. the mean chain, is always available as the alpha = 1 answer).  Passing RiskIndex's
    // geometric grid isolates the selection rule (surrogate vs. exact CVaR) from threshold placement.
    SeqGreedyCache(const GRRSketch& s, uint32_t Kmax, uint32_t Hgrid, const std::vector<double>* grid = nullptr) : sk(s), covered(s.num_sets(), 0) {
        Timer T;
        ChainResult mc = build_chain(sk, TAU_INF, Kmax, false);
        mean_chain = mc.seeds;
        sk.profile(mean_chain, X, covered, touched);
        double Gamma = *std::max_element(X.begin(), X.end()); if (Gamma <= 0) Gamma = sk.scale();
        if (grid) { for (double t : *grid) if (t != TAU_INF) taus.push_back(t); } else taus = uniform_grid(Gamma, Hgrid);
        for (double tau : taus) chains.push_back(build_chain(sk, tau, Kmax, false));
        build_sec = T.sec();
    }
    std::vector<uint32_t> query(uint32_t k, double alpha, bool select_cvar = false) {
        if (alpha >= 1.0) return std::vector<uint32_t>(mean_chain.begin(), mean_chain.begin() + k);
        double best = -std::numeric_limits<double>::infinity(); size_t bh = 0;
        for (size_t h = 0; h < chains.size(); ++h) {
            std::vector<uint32_t> S(chains[h].seeds.begin(), chains[h].seeds.begin() + std::min<size_t>(k, chains[h].seeds.size()));
            sk.profile(S, X, covered, touched);
            double v = select_cvar ? cvar_unsorted(X, alpha) : surrogate_value(X, taus[h], alpha);
            if (v > best) { best = v; bh = h; }
        }
        return std::vector<uint32_t>(chains[bh].seeds.begin(), chains[bh].seeds.begin() + std::min<size_t>(k, chains[bh].seeds.size()));
    }
};

}  // namespace ri
