// RiskIndex: Grouped Reverse-Reachability (GRR) sketch.
//
// M live-edge worlds L_1..L_M; within world m, b uniformly sampled targets
// v_{m1..mb}; R_{mj} = RR(v_{mj}; L_m).  All b RR sets of group m share the
// same world.  Setting grouped=false produces an ordinary independent-RR sketch
// of identical shape (M*b sets, one fresh world per set) that is post-hoc
// arranged into M pseudo-groups of size b -- the matched-budget baseline.
#pragma once
#include "../diffusion.hpp"

namespace ri {

struct GRRSketch {
    uint32_t n = 0, M = 0, b = 0;
    bool grouped = true;
    uint64_t L = 0;                 // total incidences sum_{m,j} |R_mj|
    std::vector<uint64_t> inv_off;  // node -> list of RR ids (global id = m*b+j)
    std::vector<uint32_t> inv_rr;
    double build_sec = 0;
    uint64_t edges_examined = 0;

    inline uint32_t num_sets() const { return M * b; }
    inline uint32_t group_of(uint32_t rrid) const { return rrid / b; }
    inline double scale() const { return (double)n / (double)b; }  // n/b
    inline uint64_t bytes() const { return inv_off.size() * 8 + inv_rr.size() * 4; }
    inline uint32_t deg(uint32_t u) const { return (uint32_t)(inv_off[u + 1] - inv_off[u]); }

    static GRRSketch build(const Diffusion& d, uint32_t M, uint32_t b, uint64_t seed, bool grouped) {
        Timer T;
        GRRSketch s; s.n = d.g.n; s.M = M; s.b = b; s.grouped = grouped;
        uint32_t n = d.g.n;
        std::vector<uint32_t> stamp(n, 0), out, queue;
        std::vector<uint32_t> nodes;          // flat node list of all RR sets
        std::vector<uint32_t> cnt(n + 1, 0);  // per node incidence count
        std::vector<uint64_t> set_off; set_off.reserve((size_t)M * b + 1); set_off.push_back(0);
        Rng target_rng(seed ^ 0x7777ULL);
        uint32_t cur = 0;
        for (uint32_t m = 0; m < M; ++m) {
            uint64_t w_group = mix64(seed) ^ mix64(m + 1);
            for (uint32_t j = 0; j < b; ++j) {
                uint32_t rrid = m * b + j;
                uint64_t w = grouped ? w_group : (mix64(seed) ^ mix64((uint64_t)rrid + 0x100000ULL));
                uint32_t t = target_rng.below(n);
                ++cur;
                d.rr_set(w, t, stamp, cur, out, queue);
                for (uint32_t u : out) { nodes.push_back(u); cnt[u + 1]++; }
                set_off.push_back(nodes.size());
            }
        }
        s.L = nodes.size();
        for (uint32_t i = 0; i < n; ++i) cnt[i + 1] += cnt[i];
        s.inv_off.assign(cnt.begin(), cnt.end());
        s.inv_rr.resize(s.L);
        std::vector<uint64_t> pos(cnt.begin(), cnt.end() - 1);
        for (uint32_t rrid = 0; rrid < M * b; ++rrid)
            for (uint64_t i = set_off[rrid]; i < set_off[rrid + 1]; ++i) s.inv_rr[pos[nodes[i]]++] = rrid;
        s.build_sec = T.sec();
        return s;
    }

    // Scenario profile of a seed set: hat X_m(S) = (n/b) * #{j : S ∩ R_mj ≠ ∅}.
    // `covered` is a scratch bitmap over RR ids (reset after use).
    void profile(const std::vector<uint32_t>& S, std::vector<double>& X, std::vector<uint8_t>& covered,
                 std::vector<uint32_t>& touched) const {
        X.assign(M, 0.0); touched.clear();
        for (uint32_t u : S)
            for (uint64_t i = inv_off[u]; i < inv_off[u + 1]; ++i) {
                uint32_t r = inv_rr[i];
                if (!covered[r]) { covered[r] = 1; touched.push_back(r); X[r / b] += 1.0; }
            }
        for (uint32_t r : touched) covered[r] = 0;
        double sc = scale();
        for (double& x : X) x *= sc;
    }
};

}  // namespace ri
