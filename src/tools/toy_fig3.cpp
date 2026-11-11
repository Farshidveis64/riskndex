// Real run of RiskIndex on the 8-node running example of Fig. 3 (M = b = 4, H = 3, Kmax = 5).
// Uses the library code unchanged (GRRSketch::build, RiskIndex::build, query); this driver only
// records what the figure shows.  Usage:
//   toy_fig3 search model=IC p=0.5 seeds=0-3000      candidate seeds that illustrate the method
//   toy_fig3 dump   model=IC p=0.5 seed=S vm=0|1     every number drawn in Fig. 3 / Fig. 2
#include "../graph.hpp"
#include "../diffusion.hpp"
#include "../grr/grr.hpp"
#include "../risk_basis/risk_basis.hpp"
#include "../index/index.hpp"
#include <cstdio>
#include <string>
#include <map>
using namespace ri;
static std::map<std::string, std::string> kv;
static std::string arg(const std::string& k, const std::string& d) { auto it = kv.find(k); return it == kv.end() ? d : it->second; }
static Graph toy_graph() {   // E8 of the figure, v1..v8 -> 0..7
    int E[12][2] = {{1,2},{1,3},{3,6},{3,4},{2,4},{4,5},{4,8},{4,7},{6,7},{8,7},{2,5},{3,7}};
    std::vector<std::pair<uint32_t,uint32_t>> e; for (auto& x : E) e.push_back({(uint32_t)x[0]-1,(uint32_t)x[1]-1});
    return Graph::from_edges(8, e);
}
static double mean_of_P(const std::vector<double>& P) { return P.back() / (P.size() - 1); }
int main(int argc, char** argv) {
    std::string cmd = argv[1];
    for (int i = 2; i < argc; ++i) { std::string s = argv[i]; auto q = s.find('='); kv[s.substr(0, q)] = s.substr(q + 1); }
    Graph g = toy_graph();
    Diffusion d(g, parse_model(arg("model", "IC")), std::stod(arg("p", "0.5")));
    const uint32_t M = 4, b = 4, H = 3, K = 5; const uint32_t k = 2; const double alpha = 0.25;
    if (cmd == "search") {
        auto rg = arg("seeds", "0-2000"); long s0 = std::stol(rg.substr(0, rg.find('-'))), s1 = std::stol(rg.substr(rg.find('-') + 1));
        bool vm = arg("vm", "0") == "1";
        for (long seed = s0; seed <= s1; ++seed) {
            GRRSketch sk = GRRSketch::build(d, M, b, (uint64_t)seed, true);
            RiskIndex idx = RiskIndex::build(sk, H, K, 0, vm);
            auto a = idx.query(k, alpha);
            uint32_t hm = H - 1;
            const auto& Pa = idx.chains[a.h].prefix_sum[k - 1]; const auto& Pm = idx.chains[hm].prefix_sum[k - 1];
            double ca = cvar_from_prefix(Pa, alpha), cm = cvar_from_prefix(Pm, alpha);
            double ma = mean_of_P(Pa), mm = mean_of_P(Pm);
            bool distinct = true;
            for (uint32_t h1 = 0; h1 < H; ++h1) for (uint32_t h2 = h1 + 1; h2 < H; ++h2) {
                auto x = std::vector<uint32_t>(idx.chains[h1].seeds.begin(), idx.chains[h1].seeds.begin() + k);
                auto y = std::vector<uint32_t>(idx.chains[h2].seeds.begin(), idx.chains[h2].seeds.begin() + k);
                std::sort(x.begin(), x.end()); std::sort(y.begin(), y.end()); if (x == y) distinct = false; }
            double ra = idx.chains[a.h].rho[k - 1], rm = idx.chains[hm].rho[k - 1];
            double second = -1; for (uint32_t h = 0; h < H; ++h) if (h != a.h) second = std::max(second, cvar_from_prefix(idx.chains[h].prefix_sum[k - 1], alpha));
            if (ca - second < std::stod(arg("margin", "0"))) continue;
            bool need_distinct = arg("distinct", "1") == "1";
            if (a.h != hm && ca > cm && mm > ma && (distinct || !need_distinct) && ra >= std::stod(arg("minrho", "0")))
                printf("seed=%ld h=%u cvar=%.3g vs mean-chain %.3g | mean %.3g vs %.3g | rho %.2f %.2f | taus %.3g %.3g\n", seed, a.h, ca, cm, ma, mm, ra, rm, idx.taus[0], idx.taus[1]);
        }
        return 0;
    }
    // dump
    uint64_t seed = std::stoull(arg("seed", "0")); bool vm = arg("vm", "0") == "1";
    GRRSketch sk = GRRSketch::build(d, M, b, seed, true);
    RiskIndex idx = RiskIndex::build(sk, H, K, 0, vm);
    printf("{\n \"model\":\"%s\",\"p\":%g,\"seed\":%llu,\"vm\":%d,\"M\":%u,\"b\":%u,\"H\":%u,\"Kmax\":%u,\n", arg("model","IC").c_str(), d.p, (unsigned long long)seed, (int)vm, M, b, H, K);
    // worlds: live edges and RR sets (replicates the sampling order of GRRSketch::build)
    printf(" \"worlds\":[\n");
    Rng target_rng(seed ^ 0x7777ULL); std::vector<uint32_t> stamp(g.n, 0), out, queue; uint32_t cur = 0;
    for (uint32_t m = 0; m < M; ++m) {
        uint64_t w = mix64(seed) ^ mix64(m + 1);
        printf("  {\"live\":[");
        bool first = true;
        for (uint32_t u = 0; u < g.n; ++u) for (uint64_t i = g.fwd_off[u]; i < g.fwd_off[u + 1]; ++i) {
            uint32_t v = g.fwd_dst[i]; bool lv = d.model == Model::LT ? (g.rev_src[g.rev_off[v] + d.lt_choice(w, v)] == u) : d.live(w, u, v);
            if (lv) { printf("%s[%u,%u]", first ? "" : ",", u + 1, v + 1); first = false; } }
        printf("],\"rr\":[");
        for (uint32_t j = 0; j < b; ++j) {
            uint32_t t = target_rng.below(g.n); ++cur; d.rr_set(w, t, stamp, cur, out, queue);
            printf("%s{\"target\":%u,\"nodes\":[", j ? "," : "", t + 1);
            for (size_t q = 0; q < out.size(); ++q) printf("%s%u", q ? "," : "", out[q] + 1);
            printf("]}");
        }
        printf("]}%s\n", m + 1 < M ? "," : "");
    }
    printf(" ],\n \"taus\":[");
    for (uint32_t h = 0; h < H; ++h) printf("%s%s", h ? "," : "", std::isinf(idx.taus[h]) ? "\"inf\"" : std::to_string(idx.taus[h]).c_str());
    printf("],\n \"chains\":[\n");
    for (uint32_t h = 0; h < H; ++h) {
        printf("  {\"seeds\":["); for (size_t i = 0; i < idx.chains[h].seeds.size(); ++i) printf("%s%u", i ? "," : "", idx.chains[h].seeds[i] + 1);
        printf("],\"cvar\":["); for (size_t kk = 0; kk < idx.chains[h].prefix_sum.size(); ++kk) printf("%s%.4f", kk ? "," : "", cvar_from_prefix(idx.chains[h].prefix_sum[kk], alpha));
        printf("],\"mean\":["); for (size_t kk = 0; kk < idx.chains[h].prefix_sum.size(); ++kk) printf("%s%.4f", kk ? "," : "", mean_of_P(idx.chains[h].prefix_sum[kk]));
        printf("],\"prefix\":["); for (size_t kk = 0; kk < idx.chains[h].prefix_sum.size(); ++kk) { printf("%s[", kk ? "," : ""); auto& P = idx.chains[h].prefix_sum[kk]; for (size_t r = 0; r < P.size(); ++r) printf("%s%.4f", r ? "," : "", P[r]); printf("]"); }
        printf("]}%s\n", h + 1 < H ? "," : "");
    }
    auto a = idx.query(k, alpha);
    // per-world profiles (unsorted, in world order) and rho of the answer and of the mean-chain prefix
    std::vector<uint8_t> cov(sk.num_sets(), 0); std::vector<uint32_t> touched; std::vector<double> X;
    auto prof = [&](const std::vector<uint32_t>& S) { sk.profile(S, X, cov, touched); std::vector<double> Y = X; double rho = vm_shrink(Y, g.n, b); return std::make_pair(X, rho); };
    std::vector<uint32_t> Sa(a.seeds), Sm(idx.chains[H - 1].seeds.begin(), idx.chains[H - 1].seeds.begin() + k);
    auto pa = prof(Sa); auto pm = prof(Sm);
    printf(" ],\n \"query\":{\"k\":%u,\"alpha\":%g,\"h\":%u,\"cert\":%.4f,\"seeds\":[", k, alpha, a.h, a.score);
    for (size_t i = 0; i < a.seeds.size(); ++i) printf("%s%u", i ? "," : "", a.seeds[i] + 1);
    printf("]},\n \"answer_profile\":["); for (size_t m = 0; m < pa.first.size(); ++m) printf("%s%.4f", m ? "," : "", pa.first[m]);
    printf("],\"answer_rho\":%.4f,\n \"meanchain_profile\":[", pa.second); for (size_t m = 0; m < pm.first.size(); ++m) printf("%s%.4f", m ? "," : "", pm.first[m]);
    printf("],\"meanchain_rho\":%.4f\n}\n", pm.second);
    return 0;
}
