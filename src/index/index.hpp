// RiskIndex: materialized prefix-profile index.
//
// Offline: H threshold-specific greedy chains over the GRR sketch, each to Kmax.
// For every prefix S_h^{(k)} the sorted scenario profile is stored as prefix sums
// P^{h,k}_0..M, which is all an online (k, alpha) query reads.
#pragma once
#include "../risk_basis/risk_basis.hpp"
#include <limits>

namespace ri {

struct ThresholdGrid {
    // Data-driven geometric grid derived only from the sketch (never from the
    // evaluation data): H-1 thresholds geometric in [tau_lo, tau_hi] plus tau=inf.
    static std::vector<double> make(const GRRSketch& sk, const ChainResult& mean_chain, uint32_t H,
                                    double alpha_min = 0.02) {
        std::vector<double> taus;
        if (H == 1) { taus.push_back(TAU_INF); return taus; }
        const auto& P1 = mean_chain.prefix_sum.front();
        const auto& PK = mean_chain.prefix_sum.back();
        size_t M = P1.size() - 1;
        // alpha_min-quantile of the 1-prefix profile (lower bound on useful VaR)
        size_t r = (size_t)std::ceil(alpha_min * M); if (r == 0) r = 1;
        double tau_lo = std::max(sk.scale(), P1[r] - P1[r - 1]);
        double tau_hi = PK[M] - PK[M - 1];  // largest scenario of the Kmax prefix
        if (tau_hi <= tau_lo) tau_hi = tau_lo * 2;
        for (uint32_t h = 0; h + 1 < H; ++h) {
            double t = (H == 2) ? std::sqrt(tau_lo * tau_hi)
                                : tau_lo * std::pow(tau_hi / tau_lo, (double)h / (double)(H - 2));
            taus.push_back(t);
        }
        taus.push_back(TAU_INF);
        return taus;
    }
};

struct RiskIndex {
    uint32_t H = 0, Kmax = 0, M = 0;   // M = number of worlds whose profile is stored
    uint32_t M_build = 0;              // worlds used by greedy (cross-fitting if < sketch M)
    std::vector<double> taus;
    std::vector<ChainResult> chains;   // chains[h].seeds, chains[h].prefix_sum[k-1]
    double build_sec = 0;
    uint64_t gain_evals = 0;

    // Exact size of the served structure per prefix: (M+1) doubles of prefix sums,
    // (M+1) doubles of prefix sums of squares when stored, rho and the mean (2
    // doubles), and one 32-bit seed id.
    bool store_sq = false;
    uint64_t bytes() const { return (uint64_t)H * Kmax * ((M + 1) * 8 * (store_sq ? 2 : 1) + 16 + 4); }

    // Algorithm 1 (chain construction + materialization) given a built sketch.
    bool vm = true;                    // variance-matched profiles
    // Selection rule (Algorithm 2).  z = 0: plain argmax of the certificate.
    // z > 0: noise margin -- a chain other than the mean chain is returned only if
    // its certificate exceeds that of the mean chain by more than z times the
    // combined tail standard error of the two certificates (requires store_sq).
    double z = 0.0;
    // mean_seq != nullptr: hybrid index -- the tau = inf chain is the given seed sequence (RIS greedy on an
    // independent RR sketch), materialized on the GRR sketch; all other chains are built on the GRR sketch.
    bool hybrid = false;
    static RiskIndex build(const GRRSketch& sk, uint32_t H, uint32_t Kmax, uint32_t M_build = 0, bool vm = true, const std::vector<double>* fixed_taus = nullptr,
                           bool store_sq = false, double z = 0.0, const std::vector<uint32_t>* mean_seq = nullptr) {
        Timer T; RiskIndex idx; idx.H = H; idx.Kmax = Kmax; idx.vm = vm; idx.store_sq = store_sq; idx.z = z;
        idx.M_build = (M_build == 0 || M_build >= sk.M) ? sk.M : M_build;
        idx.M = (idx.M_build == sk.M) ? sk.M : sk.M - idx.M_build;
        idx.hybrid = (mean_seq != nullptr);
        ChainResult mean_chain = idx.hybrid ? materialize_sequence(sk, *mean_seq, vm, store_sq)
                                            : build_chain(sk, TAU_INF, Kmax, true, idx.M_build, vm, store_sq);
        idx.gain_evals += mean_chain.gain_evals;
        idx.taus = fixed_taus ? *fixed_taus : ThresholdGrid::make(sk, mean_chain, H);
        idx.chains.resize(H);
        for (uint32_t h = 0; h < H; ++h) {
            if (idx.taus[h] == TAU_INF) idx.chains[h] = mean_chain;
            else { idx.chains[h] = build_chain(sk, idx.taus[h], Kmax, true, idx.M_build, vm, store_sq); idx.gain_evals += idx.chains[h].gain_evals; }
        }
        idx.n_nodes = sk.n; idx.b_sets = sk.b;
        idx.build_sec = T.sec();
        return idx;
    }
    // Validation: re-materialize every stored chain on an independent GRR sketch (same n, M, b), so that each
    // stored profile -- and hence each certificate -- is evaluated on worlds that did not select the chain.
    bool validated = false;
    void rematerialize(const GRRSketch& V) {
        Timer T;
        for (auto& c : chains) { double tau = c.tau; uint64_t ge = c.gain_evals; c = materialize_sequence(V, c.seeds, vm, store_sq); c.tau = tau; c.gain_evals = ge; }
        M = V.M; M_build = V.M; n_nodes = V.n; b_sets = V.b; validated = true; build_sec += T.sec();
    }
    inline uint32_t mean_h() const { for (uint32_t h = 0; h < H; ++h) if (taus[h] == TAU_INF) return h; return H - 1; }

    // Algorithm 2: online query. O(H) profile reads + O(k) to emit the seeds.
    // Answer fields (all O(1) reads of the returned entry):
    //   score : certificate = empirical lower-tail CVaR of the stored (variance-matched) profile
    //   rho   : tail-resolution ratio of that profile (diagnostic, not an error bound)
    //   raw   : raw (uncorrected) certificate, recovered exactly from (rho, mean, score) because the
    //           correction is affine: C_raw = mean - (mean - C_vm)/rho.  Defined only for rho > 0;
    //           when rho == 0 the corrected profile is constant and the raw tail cannot be recovered,
    //           so raw = NaN and raw_valid = false (the caller must not report it as a certificate).
    //   se    : tail standard error of the certificate (0 if prefix sums of squares are not stored)
    //   guarded : true if the rho-guard replaced the selected chain by the mean chain (rho_min > 0)
    struct Answer { uint32_t h; double score; double rho; double se; double se_w = 0; double raw; bool raw_valid; bool guarded; std::vector<uint32_t> seeds; };
    double n_nodes = 0, b_sets = 0;  // for the within-world standard error
    // rho-guard: if the resolution ratio of the selected (non-mean) entry is below rho_min, the sketch
    // does not resolve the tail it would be selecting on, and the mean chain is returned instead.
    double rho_min = 0.0;
    // margin_within: the z-margin uses the standard error that includes the within-world noise
    // of the tail scenarios (tail_se_within) instead of the between-world one only.
    bool margin_within = false;
    inline double se_of(uint32_t h, uint32_t k, double alpha) const {
        const auto& P = chains[h].prefix_sum[k - 1]; const auto& Q = chains[h].prefix_sq[k - 1];
        return margin_within ? tail_se_within(P, Q, alpha, chains[h].mean[k - 1], vm ? chains[h].rho[k - 1] : 1.0, n_nodes, b_sets)
                             : tail_se_from_prefix(P, Q, alpha);
    }
    inline double raw_cert(uint32_t h, uint32_t k, double c, bool* valid = nullptr) const {
        double rho = chains[h].rho[k - 1], mu = chains[h].mean[k - 1];
        if (valid) *valid = true;
        if (!vm || rho >= 1.0) return c;
        if (rho <= 0.0) { if (valid) *valid = false; return std::numeric_limits<double>::quiet_NaN(); }
        return mu - (mu - c) / rho;
    }
    // Selection shared by query(), query_score() and budget_query(): argmax of the stored certificate over the
    // chains, then the optional noise margin and the rho-guard.  Returns (h, certificate, guarded).
    struct Sel { uint32_t h; double score; bool guarded; };
    inline Sel select(uint32_t k, double alpha) const {
        Sel a; a.h = 0; a.score = -1; a.guarded = false;
        for (uint32_t h = 0; h < H; ++h) {
            if (chains[h].prefix_sum.size() < k) continue;
            double c = cvar_from_prefix(chains[h].prefix_sum[k - 1], alpha);
            if (c > a.score) { a.score = c; a.h = h; }
        }
        uint32_t hm = mean_h();
        if (z > 0 && store_sq && a.h != hm) {
            double cm = cvar_from_prefix(chains[hm].prefix_sum[k - 1], alpha);
            double s1 = se_of(a.h, k, alpha), s2 = se_of(hm, k, alpha);
            if (a.score - cm <= z * std::sqrt(s1 * s1 + s2 * s2)) { a.h = hm; a.score = cm; }
        }
        if (rho_min > 0 && a.h != hm && chains[a.h].rho[k - 1] < rho_min) {
            a.h = hm; a.score = cvar_from_prefix(chains[hm].prefix_sum[k - 1], alpha); a.guarded = true;
        }
        return a;
    }
    Answer query(uint32_t k, double alpha) const {
        Sel s = select(k, alpha);
        Answer a; a.h = s.h; a.score = s.score; a.se = 0; a.guarded = s.guarded;
        a.rho = chains[a.h].rho[k - 1];
        if (store_sq) {
            a.se = tail_se_from_prefix(chains[a.h].prefix_sum[k - 1], chains[a.h].prefix_sq[k - 1], alpha);
            a.se_w = tail_se_within(chains[a.h].prefix_sum[k - 1], chains[a.h].prefix_sq[k - 1], alpha, chains[a.h].mean[k - 1], vm ? chains[a.h].rho[k - 1] : 1.0, n_nodes, b_sets);
        }
        a.raw = raw_cert(a.h, k, a.score, &a.raw_valid);
        a.seeds.assign(chains[a.h].seeds.begin(), chains[a.h].seeds.begin() + k);
        return a;
    }
    // Budget query: smallest k <= Kmax whose returned certificate (same selection as query(), guard included)
    // reaches the target gamma.  Linear scan in k, O(H Kmax): variance-matched certificates need not be monotone
    // in k (rho_hat changes with k), so a binary search would not be exact.  Returns (k, h) or (0, H) if infeasible.
    std::pair<uint32_t, uint32_t> budget_query(double alpha, double gamma) const {
        uint32_t K = 0; for (uint32_t h = 0; h < H; ++h) K = std::max<uint32_t>(K, (uint32_t)chains[h].prefix_sum.size());
        for (uint32_t k = 1; k <= K; ++k) { Sel s = select(k, alpha); if (s.score >= gamma) return {k, s.h}; }
        return {0, H};
    }
    // Chain and certificate only (no seed copy, no standard error): identical selection to query().
    inline std::pair<uint32_t, double> query_score(uint32_t k, double alpha) const {
        Sel s = select(k, alpha); return {s.h, s.score};
    }
};

// Baseline "no materialization": chains are kept but profiles are recomputed from
// the sketch's inverted lists at query time.
struct RiskIndexNoMat {
    const GRRSketch& sk;
    const RiskIndex& idx;
    std::vector<uint8_t> covered; std::vector<uint32_t> touched; std::vector<double> X;
    RiskIndexNoMat(const GRRSketch& s, const RiskIndex& i) : sk(s), idx(i), covered(s.num_sets(), 0) {}
    std::pair<uint32_t, double> query(uint32_t k, double alpha) {
        uint32_t best = 0; double bs = -1;
        for (uint32_t h = 0; h < idx.H; ++h) {
            std::vector<uint32_t> S(idx.chains[h].seeds.begin(), idx.chains[h].seeds.begin() + k);
            sk.profile(S, X, covered, touched);
            if (idx.vm) vm_shrink(X, (double)sk.n, (double)sk.b);
            std::sort(X.begin(), X.end());
            double c = cvar_sorted(X, alpha);
            if (c > bs) { bs = c; best = h; }
        }
        return {best, bs};
    }
};

}  // namespace ri
