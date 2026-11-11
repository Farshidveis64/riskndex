// Forward layer for under-resolved small-budget tails (Step 3b of the paper).
//
// A GRR profile resolves spreads only in units of n/b.  For budgets k <= kf whose RIS-prefix tail is under-resolved
// (need gate: tail-resolution ratio R_tail(k) = b * C_0.05 / n < theta) and cheap to simulate (cost gate: estimated
// RIS-prefix spread <= sigma_max), the spreads of the candidate seed sets are small, so forward simulation of live-edge worlds
// gives their exact per-world spreads at a cost proportional to the spreads, independent of the target
// resolution n/b that limits a GRR sketch on large graphs.  The layer
//   1. screens the Pscreen best-covered nodes (RR coverage) by their exact single-seed mean and CVaR on Wsel
//      selection worlds and keeps the top-P by each statistic (plus the first kf seeds of every existing chain),
//   2. grows forward chains on the selection worlds over that pool: lazy-free greedy on the exact mean, on the exact
//      empirical CVaR_alpha (alpha in `alphas`) and on the exact truncated spreads F_tau (T thresholds, geometric
//      between the 5% quantile of the RIS 1-prefix and the 20% quantile of the RIS kf-prefix),
//   3. materializes the sorted profile of every k-prefix (k <= kf) of every forward chain and of every existing
//      chain (GRR threshold chains and the RIS chain) on Wval independent validation worlds.
// A query (k <= kf, alpha) returns the prefix with the largest exact empirical CVaR on the validation worlds.
// Because the validation worlds are independent of everything that produced the candidates, every certificate is
// an out-of-sample estimate with between-world error only (no within-world term).
#pragma once
#include "../eval.hpp"
#include "../risk_basis/risk_basis.hpp"
#include <set>
#include <map>

namespace ri {

struct ForwardLayer {
    struct Params {
        uint32_t kf = 5, Wsel = 1000, Wval = 1000, P = 60, Pscreen = 2000, T = 7;
        double sigma_max = 1000;   // cost cap: only budgets whose estimated RIS spread is at most sigma_max nodes
        double theta = 1e300;      // need gate: only budgets whose tail-resolution ratio R_tail(k) is below theta
        std::vector<double> alphas{0.05, 0.1, 0.2};
        uint64_t sel_stream = 1, val_stream = 2;
    };
    struct Cand { std::vector<uint32_t> seeds; int src; std::vector<std::vector<double>> prefix, prefix_sq; };
    struct Answer { int c = -1; double score = 0; double se = 0; std::vector<uint32_t> seeds; };
    Params prm; uint32_t W = 0, kf = 0; int ris_cand = -1;   // kf = number of budgets the layer serves (0: none)
    std::vector<Cand> cands; std::vector<double> taus;
    double build_sec = 0, screen_sec = 0, greedy_sec = 0, val_sec = 0; size_t pool_size = 0;

    // est_sigma[k-1]: estimated expected spread of the RIS k-prefix (from the RR sketch).  The layer serves the budgets
    // k <= kf with est_sigma[k-1] <= sigma_max, because the cost of forward simulation grows with the spreads.
    // rtail[k-1]: tail-resolution ratio of the RIS k-prefix, b * C_hat_0.05 / n, the number of GRR
    // resolution units n/b in the mean of its worst 5% scenarios on the validation sketch; the layer serves only the
    // leading budgets with rtail < theta (need gate) and est_sigma <= sigma_max (cost gate).  Empty rtail: no need gate.
    static ForwardLayer build(const Diffusion& d, const std::vector<uint32_t>& order, const std::vector<std::vector<uint32_t>>& chains,
                              int ris_idx, const Params& p, bool materialize, const std::vector<double>& est_sigma,
                              const std::vector<double>& rtail = {}) {
        Timer T0; ForwardLayer F; F.prm = p; F.W = p.Wval;
        uint32_t kf = 0;
        while (kf < p.kf && kf < est_sigma.size() && est_sigma[kf] <= p.sigma_max && (rtail.empty() || (kf < rtail.size() && rtail[kf] < p.theta))) ++kf;
        F.kf = kf; if (kf == 0) { F.build_sec = T0.sec(); return F; }
        const size_t A = p.alphas.size();
        // 1. screening
        Timer T1;
        uint32_t NS = std::min<uint32_t>(p.Pscreen, d.g.n);
        std::vector<std::vector<uint32_t>> singles; singles.reserve(NS);
        for (uint32_t i = 0; i < NS; ++i) singles.push_back({order[i]});
        auto ssp = evaluate_sets(d, singles, p.Wsel, p.sel_stream, materialize);
        std::vector<std::vector<double>> st(NS); std::vector<double> v(p.Wsel);
        for (uint32_t i = 0; i < NS; ++i) {
            for (uint32_t w = 0; w < p.Wsel; ++w) v[w] = ssp[i][w];
            st[i].push_back(mean_of(v)); std::sort(v.begin(), v.end());
            for (double al : p.alphas) st[i].push_back(cvar_sorted(v, al));
        }
        std::set<uint32_t> ps;
        for (size_t c = 0; c <= A; ++c) {
            std::vector<uint32_t> idx(NS); for (uint32_t i = 0; i < NS; ++i) idx[i] = i;
            uint32_t top = std::min<uint32_t>(p.P, NS);
            std::partial_sort(idx.begin(), idx.begin() + top, idx.end(), [&](uint32_t x, uint32_t y) { return st[x][c] > st[y][c]; });
            for (uint32_t i = 0; i < top; ++i) ps.insert(order[idx[i]]);
        }
        for (const auto& ch : chains) for (uint32_t j = 0; j < kf && j < ch.size(); ++j) ps.insert(ch[j]);
        std::vector<uint32_t> pool(ps.begin(), ps.end()); F.pool_size = pool.size();
        F.screen_sec = T1.sec();
        // 2. tau grid from the RIS prefixes on the selection worlds
        Timer T2;
        std::vector<std::vector<uint32_t>> rp;
        for (uint32_t k = 1; k <= kf; ++k) rp.push_back(std::vector<uint32_t>(chains[ris_idx].begin(), chains[ris_idx].begin() + k));
        auto rsp = evaluate_sets(d, rp, p.Wsel, p.sel_stream, materialize);
        auto quant = [](const std::vector<uint32_t>& x, double q) { std::vector<double> z(x.begin(), x.end()); std::sort(z.begin(), z.end()); return z[(size_t)std::floor(q * (z.size() - 1))]; };
        double tlo = std::max(1.0, quant(rsp[0], 0.05)), thi = std::max(tlo * 1.01, quant(rsp[kf - 1], 0.2));
        for (uint32_t t = 0; t < p.T; ++t) F.taus.push_back(p.T == 1 ? tlo : tlo * std::pow(thi / tlo, (double)t / (p.T - 1)));
        // 3. forward chains: objective 0 = mean, 1..A = CVaR_alpha, A+1.. = F_tau; one evaluation pass per step
        const size_t NO = 1 + A + F.taus.size();
        std::vector<std::vector<uint32_t>> cur(NO);
        for (uint32_t step = 1; step <= kf; ++step) {
            // candidate sets of all objectives; identical sets (objectives that share a prefix) are evaluated once
            std::vector<std::vector<uint32_t>> cand; std::vector<std::pair<size_t, uint32_t>> who; std::vector<size_t> slot;
            std::map<std::vector<uint32_t>, size_t> seen;
            for (size_t o = 0; o < NO; ++o) for (uint32_t u : pool) {
                if (std::find(cur[o].begin(), cur[o].end(), u) != cur[o].end()) continue;
                std::vector<uint32_t> S = cur[o]; S.push_back(u);
                std::vector<uint32_t> key = S; std::sort(key.begin(), key.end());
                auto it = seen.find(key);
                if (it == seen.end()) { seen.emplace(std::move(key), cand.size()); slot.push_back(cand.size()); cand.push_back(std::move(S)); }
                else slot.push_back(it->second);
                who.push_back({o, u});
            }
            auto sp = evaluate_sets(d, cand, p.Wsel, p.sel_stream, materialize);
            std::vector<double> best(NO, -1.0), bestm(NO, -1.0); std::vector<uint32_t> arg(NO, UINT32_MAX);
            for (size_t c = 0; c < who.size(); ++c) {
                size_t o = who[c].first;
                for (uint32_t w = 0; w < p.Wsel; ++w) v[w] = sp[slot[c]][w];
                double mu = mean_of(v), sc = mu;
                if (o >= 1 && o <= A) { std::sort(v.begin(), v.end()); sc = cvar_sorted(v, p.alphas[o - 1]); }
                else if (o > A) { double tau = F.taus[o - 1 - A]; sc = 0; for (double x : v) sc += std::min(x, tau); sc /= p.Wsel; }
                if (sc > best[o] + 1e-12 || (std::fabs(sc - best[o]) <= 1e-12 && mu > bestm[o])) { best[o] = sc; bestm[o] = mu; arg[o] = who[c].second; }
            }
            for (size_t o = 0; o < NO; ++o) cur[o].push_back(arg[o]);
        }
        F.greedy_sec = T2.sec();
        // 4. candidates: forward chains, then the existing chains (the RIS chain among them)
        for (size_t o = 0; o < NO; ++o) F.cands.push_back({cur[o], (int)o, {}, {}});
        for (size_t h = 0; h < chains.size(); ++h) {
            if ((int)h == ris_idx) F.ris_cand = (int)F.cands.size();
            F.cands.push_back({std::vector<uint32_t>(chains[h].begin(), chains[h].begin() + std::min<size_t>(kf, chains[h].size())), 100 + (int)h, {}, {}});
        }
        // 5. validation profiles of every k-prefix
        Timer T3;
        std::vector<std::vector<uint32_t>> vs; std::vector<std::pair<size_t, uint32_t>> vk;
        for (size_t c = 0; c < F.cands.size(); ++c) for (uint32_t k = 1; k <= F.cands[c].seeds.size(); ++k) {
            vs.push_back(std::vector<uint32_t>(F.cands[c].seeds.begin(), F.cands[c].seeds.begin() + k)); vk.push_back({c, k});
        }
        auto vsp = evaluate_sets(d, vs, p.Wval, p.val_stream, materialize);
        for (auto& c : F.cands) { c.prefix.resize(c.seeds.size()); c.prefix_sq.resize(c.seeds.size()); }
        std::vector<double> z(p.Wval);
        for (size_t i = 0; i < vs.size(); ++i) {
            for (uint32_t w = 0; w < p.Wval; ++w) z[w] = vsp[i][w];
            std::sort(z.begin(), z.end());
            auto& P1 = F.cands[vk[i].first].prefix[vk[i].second - 1]; auto& Q1 = F.cands[vk[i].first].prefix_sq[vk[i].second - 1];
            P1.assign(p.Wval + 1, 0.0); Q1.assign(p.Wval + 1, 0.0);
            for (uint32_t w = 0; w < p.Wval; ++w) { P1[w + 1] = P1[w] + z[w]; Q1[w + 1] = Q1[w] + z[w] * z[w]; }
        }
        F.val_sec = T3.sec(); F.build_sec = T0.sec();
        return F;
    }

    // exact empirical CVaR of the best k-prefix on the validation worlds (ties: the RIS chain, then the larger mean)
    bool serves(uint32_t k) const { return k >= 1 && k <= kf; }
    Answer query(uint32_t k, double alpha) const {
        Answer best; double bm = -1;
        for (size_t c = 0; c < cands.size(); ++c) {
            if (cands[c].seeds.size() < k) continue;
            const auto& P1 = cands[c].prefix[k - 1];
            double sc = cvar_from_prefix(P1, alpha), mu = P1.back() / W;
            bool better = sc > best.score + 1e-9 || (std::fabs(sc - best.score) <= 1e-9 && ((int)c == ris_cand || (best.c != ris_cand && mu > bm)));
            if (best.c < 0 || better) { best.c = (int)c; best.score = sc; bm = mu; }
        }
        const auto& P1 = cands[best.c].prefix[k - 1]; const auto& Q1 = cands[best.c].prefix_sq[k - 1];
        size_t r = (size_t)std::ceil(alpha * W - 1e-12); if (r < 1) r = 1; if (r > W) r = W;
        double m1 = P1[r] / r, m2 = Q1[r] / r; best.se = std::sqrt(std::max(0.0, m2 - m1 * m1) / r);
        best.seeds.assign(cands[best.c].seeds.begin(), cands[best.c].seeds.begin() + k);
        return best;
    }
    // Chain index and certificate only: identical selection to query().
    std::pair<int, double> query_score(uint32_t k, double alpha) const { Answer a = query(k, alpha); return {a.c, a.score}; }
    size_t bytes() const { size_t s = 0; for (auto& c : cands) s += c.seeds.size() * 4 + c.prefix.size() * (W + 1) * 16; return s; }
};

// Routed (adaptive) interfaces.  Routing is fixed at build time (fl.kf): a budget k <= kf is served from the forward
// directory (exact per-world spreads on W worlds, no rho-guard: the profiles carry no within-world sampling noise and
// rho_hat = 1 by convention), every other budget from the GRR directory (M worlds, variance-matched certificates,
// rho-guard).  query, query_score and budget_query share this routing and the per-directory selection.
struct AdaptiveAnswer { bool forward = false; double score = 0, se = 0, rho = 1; bool guarded = false; std::vector<uint32_t> seeds; };
inline AdaptiveAnswer adaptive_query(const RiskIndex& idx, const ForwardLayer& fl, uint32_t k, double alpha) {
    AdaptiveAnswer a;
    if (fl.serves(k)) { auto f = fl.query(k, alpha); a.forward = true; a.score = f.score; a.se = f.se; a.rho = 1.0; a.seeds = f.seeds; }
    else { auto g = idx.query(k, alpha); a.score = g.score; a.se = g.se; a.rho = g.rho; a.guarded = g.guarded; a.seeds = g.seeds; }
    return a;
}
inline double adaptive_query_score(const RiskIndex& idx, const ForwardLayer& fl, uint32_t k, double alpha) {
    return fl.serves(k) ? fl.query_score(k, alpha).second : idx.query_score(k, alpha).second;
}
// Smallest k <= Kmax whose routed certificate reaches gamma; O(Kmax * max(H, H_F)).  Returns 0 if infeasible.
inline uint32_t adaptive_budget_query(const RiskIndex& idx, const ForwardLayer& fl, double alpha, double gamma) {
    uint32_t K = 0; for (uint32_t h = 0; h < idx.H; ++h) K = std::max<uint32_t>(K, (uint32_t)idx.chains[h].prefix_sum.size());
    for (uint32_t k = 1; k <= K; ++k) if (adaptive_query_score(idx, fl, k, alpha) >= gamma) return k;
    return 0;
}

}  // namespace ri
