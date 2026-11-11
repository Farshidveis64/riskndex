// RiskIndex: directed graph in CSR form (forward + reverse) with stable edge ids.
#pragma once
#include "common.hpp"
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <iostream>
#include <cstdio>

namespace ri {

struct Graph {
    uint32_t n = 0;
    uint64_t m = 0;
    // reverse CSR: for node v, in-edges are rev_off[v]..rev_off[v+1), giving src
    std::vector<uint64_t> rev_off;
    std::vector<uint32_t> rev_src;
    // forward CSR: for node u, out-edges fwd_off[u]..fwd_off[u+1), giving dst
    std::vector<uint64_t> fwd_off;
    std::vector<uint32_t> fwd_dst;
    // An edge (u,v) is identified by the 64-bit key (u<<32)|v, used to derive its
    // scenario-consistent state in a world; no explicit edge-id arrays are stored.
    static inline uint64_t ekey(uint32_t u, uint32_t v) { return ((uint64_t)u << 32) | v; }

    inline uint32_t indeg(uint32_t v) const { return (uint32_t)(rev_off[v + 1] - rev_off[v]); }
    inline uint32_t outdeg(uint32_t u) const { return (uint32_t)(fwd_off[u + 1] - fwd_off[u]); }

    // Build from edge list (u->v).  Duplicates and self loops removed.
    static Graph from_edges(uint32_t n, std::vector<std::pair<uint32_t, uint32_t>>& edges) {
        std::sort(edges.begin(), edges.end());
        edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
        std::vector<std::pair<uint32_t, uint32_t>> clean;
        clean.reserve(edges.size());
        for (auto& e : edges) if (e.first != e.second && e.first < n && e.second < n) clean.push_back(e);
        edges.swap(clean);
        Graph g; g.n = n; g.m = edges.size();
        // edge id = position in (src-sorted) list -> forward CSR is trivial
        g.fwd_off.assign(n + 1, 0);
        for (auto& e : edges) g.fwd_off[e.first + 1]++;
        for (uint32_t i = 0; i < n; ++i) g.fwd_off[i + 1] += g.fwd_off[i];
        g.fwd_dst.resize(g.m);
        for (uint64_t i = 0; i < g.m; ++i) g.fwd_dst[i] = edges[i].second;
        // reverse CSR
        g.rev_off.assign(n + 1, 0);
        for (auto& e : edges) g.rev_off[e.second + 1]++;
        for (uint32_t i = 0; i < n; ++i) g.rev_off[i + 1] += g.rev_off[i];
        g.rev_src.resize(g.m);
        std::vector<uint64_t> pos(g.rev_off.begin(), g.rev_off.end() - 1);
        for (uint64_t i = 0; i < g.m; ++i) {
            uint32_t v = edges[i].second; uint64_t p = pos[v]++;
            g.rev_src[p] = edges[i].first;
        }
        std::vector<std::pair<uint32_t, uint32_t>>().swap(edges);
        return g;
    }

    // Load a SNAP-style edge list ("# comments", "u v" per line).  If undirected,
    // every edge is inserted in both directions.  Node ids are remapped to 0..n-1
    // in order of first appearance; isolated ids never appear in the file.
    static Graph load_snap(const std::string& path, bool undirected) {
        std::ifstream in(path);
        if (!in) throw std::runtime_error("cannot open " + path);
        std::unordered_map<uint64_t, uint32_t> id;
        id.reserve(1 << 20);
        std::vector<std::pair<uint32_t, uint32_t>> edges;
        std::string line;
        auto get = [&](uint64_t x) {
            auto it = id.find(x);
            if (it != id.end()) return it->second;
            uint32_t k = (uint32_t)id.size(); id.emplace(x, k); return k;
        };
        while (std::getline(in, line)) {
            if (line.empty() || line[0] == '#' || line[0] == '%') continue;
            uint64_t a, b;
            if (sscanf(line.c_str(), "%lu %lu", &a, &b) != 2) continue;
            uint32_t u = get(a), v = get(b);
            if (u == v) continue;
            edges.emplace_back(u, v);
            if (undirected) edges.emplace_back(v, u);
        }
        return from_edges((uint32_t)id.size(), edges);
    }

    // ---- Size sweep: padding with isolated nodes (size sweep at fixed influence structure) ----
    // The influence structure of `base` is unchanged; n grows to factor * n_base.  A uniform target of a GRR sketch
    // now falls into the region that a seed set can reach with probability n_base / n, as on a large graph whose
    // seed sets reach only a small part of it.
    static Graph pad(const Graph& base, uint32_t factor) {
        Graph g = base; uint64_t n = (uint64_t)base.n * factor;
        if (n >= 0xFFFFFFFFULL) throw std::runtime_error("pad: too many nodes");
        g.n = (uint32_t)n;
        g.fwd_off.resize(n + 1, base.fwd_off.back()); g.rev_off.resize(n + 1, base.rev_off.back());
        return g;
    }

    // ---- Size sweep: replication of a real graph (size sweep at fixed local structure) ----
    // `copies` disjoint copies of `base` plus, per copy, round(inter_frac * m_base) random inter-copy edges (each
    // inserted in both directions if `undirected`).  Copy c maps node v to c * n_base + v.
    static Graph replicate(const Graph& base, uint32_t copies, double inter_frac, bool undirected, uint64_t seed) {
        uint64_t n0 = base.n; uint64_t n = n0 * copies;
        if (n >= 0xFFFFFFFFULL) throw std::runtime_error("replicate: too many nodes");
        std::vector<std::pair<uint32_t, uint32_t>> edges;
        uint64_t per = (uint64_t)std::llround(inter_frac * (double)base.m);
        edges.reserve(base.m * copies + per * copies * (undirected ? 2 : 1));
        for (uint32_t c = 0; c < copies; ++c)
            for (uint32_t u = 0; u < base.n; ++u)
                for (uint64_t i = base.fwd_off[u]; i < base.fwd_off[u + 1]; ++i)
                    edges.emplace_back((uint32_t)(c * n0 + u), (uint32_t)(c * n0 + base.fwd_dst[i]));
        Rng r(seed ^ 0x5EEDULL);
        if (copies > 1)
            for (uint32_t c = 0; c < copies; ++c)
                for (uint64_t e = 0; e < per; ++e) {
                    uint32_t c2 = r.below(copies - 1); if (c2 >= c) ++c2;
                    uint32_t a = (uint32_t)(c * n0 + r.below(base.n)), b = (uint32_t)(c2 * n0 + r.below(base.n));
                    edges.emplace_back(a, b); if (undirected) edges.emplace_back(b, a);
                }
        return from_edges((uint32_t)n, edges);
    }

    // ---- synthetic generators (undirected skeleton -> both directions) -------
    static Graph erdos_renyi(uint32_t n, double avg_deg, uint64_t seed) {
        Rng rng(seed); double p = avg_deg / (double)(n - 1);
        std::vector<std::pair<uint32_t, uint32_t>> e;
        // geometric skipping for sparse ER
        uint64_t total = (uint64_t)n * (n - 1) / 2, idx = 0;
        double lq = std::log(1.0 - p);
        while (true) {
            double u = rng.uniform(); if (u <= 0) u = 1e-300;
            uint64_t skip = (uint64_t)std::floor(std::log(u) / lq);
            idx += skip; if (idx >= total) break;
            // map idx -> (i,j) with i<j
            uint64_t i = (uint64_t)((1 + std::sqrt(1.0 + 8.0 * idx)) / 2.0);
            while (i * (i - 1) / 2 > idx) --i;
            while ((i + 1) * i / 2 <= idx) ++i;
            uint64_t j = idx - i * (i - 1) / 2;
            e.emplace_back((uint32_t)j, (uint32_t)i); e.emplace_back((uint32_t)i, (uint32_t)j);
            idx += 1;
        }
        return from_edges(n, e);
    }
    static Graph barabasi_albert(uint32_t n, uint32_t mm, uint64_t seed) {
        Rng rng(seed);
        std::vector<uint32_t> targets; targets.reserve((size_t)2 * n * mm);
        std::vector<std::pair<uint32_t, uint32_t>> e;
        for (uint32_t i = 0; i < mm; ++i) targets.push_back(i);
        for (uint32_t v = mm; v < n; ++v) {
            std::vector<uint32_t> chosen;
            while (chosen.size() < mm) {
                uint32_t t = targets[rng.below((uint32_t)targets.size())];
                if (std::find(chosen.begin(), chosen.end(), t) == chosen.end()) chosen.push_back(t);
            }
            for (uint32_t t : chosen) { e.emplace_back(v, t); e.emplace_back(t, v); targets.push_back(t); targets.push_back(v); }
        }
        return from_edges(n, e);
    }
    static Graph watts_strogatz(uint32_t n, uint32_t k, double beta, uint64_t seed) {
        Rng rng(seed);
        std::vector<std::pair<uint32_t, uint32_t>> e;
        for (uint32_t i = 0; i < n; ++i)
            for (uint32_t d = 1; d <= k / 2; ++d) {
                uint32_t j = (i + d) % n;
                if (rng.uniform() < beta) { do { j = rng.below(n); } while (j == i); }
                e.emplace_back(i, j); e.emplace_back(j, i);
            }
        return from_edges(n, e);
    }
    static Graph sbm(uint32_t n, uint32_t blocks, double p_in, double p_out, uint64_t seed) {
        Rng rng(seed);
        std::vector<std::pair<uint32_t, uint32_t>> e;
        for (uint32_t i = 0; i < n; ++i)
            for (uint32_t j = i + 1; j < n; ++j) {
                double p = (i % blocks == j % blocks) ? p_in : p_out;
                if (rng.uniform() < p) { e.emplace_back(i, j); e.emplace_back(j, i); }
            }
        return from_edges(n, e);
    }
};

}  // namespace ri
