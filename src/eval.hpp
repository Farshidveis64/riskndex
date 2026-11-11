// RiskIndex: independent evaluation on fresh live-edge worlds, and exhaustive
// oracle search on small graphs (n <= 64, bitmask reachability).
#pragma once
#include "diffusion.hpp"
#include <array>

namespace ri {

// Realized spreads of several seed sets on R fresh worlds (stream = eval_seed).
// Returns spreads[set][world].  Worlds are shared across sets (common random numbers).
inline std::vector<std::vector<uint32_t>> evaluate_sets(const Diffusion& d, const std::vector<std::vector<uint32_t>>& sets,
                                                        uint32_t R, uint64_t eval_seed, bool materialize) {
    const Graph& g = d.g; uint32_t n = g.n;
    std::vector<std::vector<uint32_t>> out(sets.size(), std::vector<uint32_t>(R, 0));
    std::vector<uint32_t> stamp(n, 0), queue; uint32_t cur = 0;
    for (uint32_t w = 0; w < R; ++w) {
        uint64_t ws = mix64(eval_seed) ^ mix64(0xE0E0E0E0ULL + w);
        if (materialize) {
            LiveWorld L(d, ws);
            for (size_t s = 0; s < sets.size(); ++s) { ++cur; out[s][w] = L.spread(sets[s], stamp, cur, queue); }
        } else {
            for (size_t s = 0; s < sets.size(); ++s) {
                ++cur; queue.clear();
                for (uint32_t v : sets[s]) if (stamp[v] != cur) { stamp[v] = cur; queue.push_back(v); }
                size_t head = 0;
                while (head < queue.size()) {
                    uint32_t u = queue[head++];
                    for (uint64_t i = g.fwd_off[u]; i < g.fwd_off[u + 1]; ++i) {
                        uint32_t v = g.fwd_dst[i];
                        if (stamp[v] == cur) continue;
                        bool live;
                        if (d.model == Model::LT) live = (g.rev_src[g.rev_off[v] + d.lt_choice(ws, v)] == u);
                        else live = hash_uniform(ws, Graph::ekey(u, v)) < d.edge_prob(v);
                        if (live) { stamp[v] = cur; queue.push_back(v); }
                    }
                }
                out[s][w] = (uint32_t)queue.size();
            }
        }
    }
    return out;
}

// ---- exhaustive oracle on n <= 64 -------------------------------------------
struct Oracle {
    uint32_t n, R;
    std::vector<std::vector<uint64_t>> reach;  // reach[w][v] bitmask of nodes reached from v in world w
    Oracle(const Diffusion& d, uint32_t R_, uint64_t seed) : n(d.g.n), R(R_) {
        if (n > 64) throw std::runtime_error("oracle requires n <= 64");
        reach.assign(R, std::vector<uint64_t>(n, 0));
        std::vector<uint32_t> stamp(n, 0), queue; uint32_t cur = 0;
        for (uint32_t w = 0; w < R; ++w) {
            uint64_t ws = mix64(seed) ^ mix64(0x0A0A0A0AULL + w);
            LiveWorld L(d, ws);
            for (uint32_t v = 0; v < n; ++v) {
                ++cur; std::vector<uint32_t> S{v}; L.spread(S, stamp, cur, queue);
                uint64_t mask = 0; for (uint32_t x : queue) mask |= (1ULL << x);
                reach[w][v] = mask;
            }
        }
    }
    // spreads of a seed set across the R worlds
    void spreads(const std::vector<uint32_t>& S, std::vector<double>& out) const {
        out.resize(R);
        for (uint32_t w = 0; w < R; ++w) { uint64_t m = 0; for (uint32_t v : S) m |= reach[w][v]; out[w] = (double)__builtin_popcountll(m); }
    }
    // exhaustive argmax of empirical CVaR_alpha over all k-subsets, for several alphas
    std::vector<std::vector<uint32_t>> best_sets(uint32_t k, const std::vector<double>& alphas, std::vector<double>& best_vals) const {
        std::vector<std::vector<uint32_t>> best(alphas.size()); best_vals.assign(alphas.size(), -1);
        std::vector<uint32_t> S(k); for (uint32_t i = 0; i < k; ++i) S[i] = i;
        std::vector<double> sp(R);
        std::vector<uint64_t> masks(R);
        while (true) {
            for (uint32_t w = 0; w < R; ++w) { uint64_t m = 0; for (uint32_t v : S) m |= reach[w][v]; sp[w] = (double)__builtin_popcountll(m); }
            std::sort(sp.begin(), sp.end());
            for (size_t a = 0; a < alphas.size(); ++a) {
                double c = cvar_sorted(sp, alphas[a]);
                if (c > best_vals[a]) { best_vals[a] = c; best[a] = S; }
            }
            // next combination
            int i = (int)k - 1;
            while (i >= 0 && S[i] == n - k + i) --i;
            if (i < 0) break;
            S[i]++; for (uint32_t j = i + 1; j < k; ++j) S[j] = S[j - 1] + 1;
        }
        return best;
    }
};

}  // namespace ri
