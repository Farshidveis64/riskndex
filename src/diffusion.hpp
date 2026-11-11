// RiskIndex: triggering-model diffusion with scenario-consistent live-edge worlds.
//
// A "world" is identified by a 64-bit stream id w.  The live/blocked state of an
// edge (or the triggering choice of a node under LT) is a deterministic function
// of (w, edge id), so any number of reverse traversals performed in the same
// world observe the same live-edge realization L_w without materializing it.
#pragma once
#include "graph.hpp"

namespace ri {

enum class Model { IC, WC, LT };

static inline const char* model_name(Model m) {
    switch (m) { case Model::IC: return "IC"; case Model::WC: return "WC"; default: return "LT"; }
}
static inline Model parse_model(const std::string& s) {
    if (s == "IC") return Model::IC;
    if (s == "WC") return Model::WC;
    if (s == "LT") return Model::LT;
    throw std::runtime_error("unknown model " + s);
}

struct Diffusion {
    const Graph& g;
    Model model;
    double p;  // uniform IC probability (only for Model::IC)
    Diffusion(const Graph& g_, Model m, double p_) : g(g_), model(m), p(p_) {}

    inline double edge_prob(uint32_t v) const {  // probability of an in-edge of v
        return model == Model::IC ? p : 1.0 / (double)g.indeg(v);
    }
    // IC / WC: is edge (u,v) live in world w?
    inline bool live(uint64_t w, uint32_t u, uint32_t v) const {
        return hash_uniform(w, Graph::ekey(u, v)) < edge_prob(v);
    }
    // LT (uniform weights 1/indeg): index of the single in-edge chosen by v in world w
    inline uint32_t lt_choice(uint64_t w, uint32_t v) const {
        uint32_t d = g.indeg(v);
        double r = hash_uniform(w, 0x8000000000000000ULL ^ (uint64_t)v);
        uint32_t idx = (uint32_t)(r * d); if (idx >= d) idx = d - 1;
        return idx;
    }

    // Reverse-reachable set of target t in world w.  `stamp` is a per-node visit
    // array; `cur` is the stamp value for this traversal.  Appends nodes to out.
    void rr_set(uint64_t w, uint32_t t, std::vector<uint32_t>& stamp, uint32_t cur,
                std::vector<uint32_t>& out, std::vector<uint32_t>& queue) const {
        out.clear(); queue.clear();
        stamp[t] = cur; queue.push_back(t); out.push_back(t);
        size_t head = 0;
        while (head < queue.size()) {
            uint32_t v = queue[head++];
            uint64_t a = g.rev_off[v], b = g.rev_off[v + 1];
            if (model == Model::LT) {
                if (b == a) continue;
                uint64_t idx = a + lt_choice(w, v);
                uint32_t u = g.rev_src[idx];
                if (stamp[u] != cur) { stamp[u] = cur; queue.push_back(u); out.push_back(u); }
            } else {
                double pe = edge_prob(v);
                for (uint64_t i = a; i < b; ++i) {
                    uint32_t u = g.rev_src[i];
                    if (stamp[u] == cur) continue;
                    if (hash_uniform(w, Graph::ekey(u, v)) < pe) { stamp[u] = cur; queue.push_back(u); out.push_back(u); }
                }
            }
        }
    }
};

// A materialized live-edge world used for *independent evaluation*: forward CSR
// of live edges only, so that the realized spread of many seed sets can be read
// off the same realization (common random numbers across methods).
struct LiveWorld {
    uint32_t n;
    std::vector<uint64_t> off;
    std::vector<uint32_t> dst;

    LiveWorld(const Diffusion& d, uint64_t w) {
        const Graph& g = d.g; n = g.n;
        std::vector<uint64_t> cnt(n + 1, 0);
        std::vector<uint8_t> is_live;
        if (d.model == Model::LT) {
            // each node v with indeg>0 has exactly one live in-edge
            off.assign(n + 1, 0);
            std::vector<uint32_t> src_of(n, UINT32_MAX);
            for (uint32_t v = 0; v < n; ++v) {
                if (g.indeg(v) == 0) continue;
                uint32_t u = g.rev_src[g.rev_off[v] + d.lt_choice(w, v)];
                src_of[v] = u; cnt[u + 1]++;
            }
            for (uint32_t i = 0; i < n; ++i) cnt[i + 1] += cnt[i];
            off = cnt; dst.resize(off[n]);
            std::vector<uint32_t> pos(off.begin(), off.end() - 1);
            for (uint32_t v = 0; v < n; ++v) if (src_of[v] != UINT32_MAX) dst[pos[src_of[v]]++] = v;
        } else {
            is_live.assign(g.m, 0);
            for (uint32_t u = 0; u < n; ++u) {
                for (uint64_t i = g.fwd_off[u]; i < g.fwd_off[u + 1]; ++i) {
                    uint32_t v = g.fwd_dst[i];
                    if (hash_uniform(w, Graph::ekey(u, v)) < d.edge_prob(v)) { is_live[i] = 1; cnt[u + 1]++; }
                }
            }
            for (uint32_t i = 0; i < n; ++i) cnt[i + 1] += cnt[i];
            off = cnt; dst.resize(off[n]);
            uint64_t p = 0;
            for (uint32_t u = 0; u < n; ++u)
                for (uint64_t i = g.fwd_off[u]; i < g.fwd_off[u + 1]; ++i) if (is_live[i]) dst[p++] = g.fwd_dst[i];
        }
    }
    // realized spread |Reach_L(S)|
    uint32_t spread(const std::vector<uint32_t>& S, std::vector<uint32_t>& stamp, uint32_t cur,
                    std::vector<uint32_t>& queue) const {
        queue.clear();
        for (uint32_t s : S) if (stamp[s] != cur) { stamp[s] = cur; queue.push_back(s); }
        size_t head = 0;
        while (head < queue.size()) {
            uint32_t u = queue[head++];
            for (uint64_t i = off[u]; i < off[u + 1]; ++i) {
                uint32_t v = dst[i];
                if (stamp[v] != cur) { stamp[v] = cur; queue.push_back(v); }
            }
        }
        return (uint32_t)queue.size();
    }
};

}  // namespace ri
