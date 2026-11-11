// RiskIndex experiment driver.  All experiments are deterministic given seeds.
//
//   riskindex oracle   graph=... model=... [k=3] [alphas=...] ...
//   riskindex fidelity graph=... model=... M=256 bs=1,4,16,64,256 ...
//   riskindex quality  graph=... model=... M=256 b=2048 H=8 Kmax=100 ks=10,25,50,100 alphas=... seeds=0,1,2
//   riskindex scale    graph=... model=... (build/latency/workload/break-even)
//
#include "graph.hpp"
#include "diffusion.hpp"
#include "grr/grr.hpp"
#include "risk_basis/risk_basis.hpp"
#include "index/index.hpp"
#include "baselines/baselines.hpp"
#include "eval.hpp"
#include "forward/forward.hpp"
#include <map>
#include <memory>
#include <set>
#include <iostream>
#include <iomanip>

using namespace ri;

// ---------------------------------------------------------------- arg parsing
struct Args {
    std::map<std::string, std::string> kv;
    std::string get(const std::string& k, const std::string& def = "") const { auto it = kv.find(k); return it == kv.end() ? def : it->second; }
    double getd(const std::string& k, double def) const { auto it = kv.find(k); return it == kv.end() ? def : std::stod(it->second); }
    long geti(const std::string& k, long def) const { auto it = kv.find(k); return it == kv.end() ? def : std::stol(it->second); }
    bool has(const std::string& k) const { return kv.count(k) > 0; }
};
static std::vector<std::string> split(const std::string& s, char c) {
    std::vector<std::string> out; std::string cur;
    for (char ch : s) { if (ch == c) { out.push_back(cur); cur.clear(); } else cur += ch; }
    out.push_back(cur); return out;
}
static std::vector<double> dlist(const std::string& s) { std::vector<double> o; for (auto& t : split(s, ',')) if (!t.empty()) o.push_back(std::stod(t)); return o; }
static std::vector<long> ilist(const std::string& s) { std::vector<long> o; for (auto& t : split(s, ',')) if (!t.empty()) o.push_back(std::stol(t)); return o; }

// graph=snap:PATH:directed|undirected | ba:n:m:seed | er:n:deg:seed | ws:n:k:beta:seed | sbm:n:blocks:pin:pout:seed
static Graph load_graph(const std::string& spec, std::string& family) {
    auto f = split(spec, ':');
    family = f[0];
    if (f[0] == "snap") return Graph::load_snap(f[1], f[2] == "undirected");
    if (f[0] == "ba") return Graph::barabasi_albert((uint32_t)std::stol(f[1]), (uint32_t)std::stol(f[2]), std::stoull(f[3]));
    if (f[0] == "er") return Graph::erdos_renyi((uint32_t)std::stol(f[1]), std::stod(f[2]), std::stoull(f[3]));
    if (f[0] == "ws") return Graph::watts_strogatz((uint32_t)std::stol(f[1]), (uint32_t)std::stol(f[2]), std::stod(f[3]), std::stoull(f[4]));
    if (f[0] == "pad") { Graph b0 = Graph::load_snap(f[1], f[2] == "undirected"); return Graph::pad(b0, (uint32_t)std::stol(f[3])); }
    if (f[0] == "rep") { Graph b0 = Graph::load_snap(f[1], f[2] == "undirected"); return Graph::replicate(b0, (uint32_t)std::stol(f[3]), std::stod(f[4]), f[2] == "undirected", std::stoull(f[5])); }
    if (f[0] == "sbm") return Graph::sbm((uint32_t)std::stol(f[1]), (uint32_t)std::stol(f[2]), std::stod(f[3]), std::stod(f[4]), std::stoull(f[5]));
    throw std::runtime_error("bad graph spec " + spec);
}

static std::string vec_str(const std::vector<uint32_t>& v) {
    std::string s; for (size_t i = 0; i < v.size(); ++i) { if (i) s += ' '; s += std::to_string(v[i]); } return s;
}
static std::string vec_str(const std::vector<double>& v) {
    std::string s; char buf[32]; for (size_t i = 0; i < v.size(); ++i) { if (i) s += ' '; snprintf(buf, 32, "%g", v[i]); s += buf; } return s;
}
static double percentile(std::vector<double> v, double p) {
    if (v.empty()) return 0;
    std::sort(v.begin(), v.end());
    double idx = p * (v.size() - 1); size_t lo = (size_t)idx; double fr = idx - lo;
    return lo + 1 < v.size() ? v[lo] * (1 - fr) + v[lo + 1] * fr : v[lo];
}

// ------------------------------------------------------------- b selection rule
// Within-world resolution rule (Theorem 2, relative form): choose b so that the
// relative standard error sqrt(n/(b*sigma)) of a per-world estimate is about
// 1/sqrt(c) for a reference seed set (the top-10 out-degree nodes), using a
// 100-cascade pilot on a separate random stream.  b is rounded up to a power of
// two and clamped to [64, bmax].
static uint32_t auto_b(const Diffusion& d, double c, uint32_t bmax) {
    const Graph& g = d.g;
    std::vector<uint32_t> by_deg(g.n); for (uint32_t v = 0; v < g.n; ++v) by_deg[v] = v;
    std::partial_sort(by_deg.begin(), by_deg.begin() + std::min<uint32_t>(10, g.n), by_deg.end(), [&](uint32_t x, uint32_t y) { return g.outdeg(x) > g.outdeg(y); });
    std::vector<std::vector<uint32_t>> sets{std::vector<uint32_t>(by_deg.begin(), by_deg.begin() + std::min<uint32_t>(10, g.n))};
    auto sp = evaluate_sets(d, sets, 100, 123456789, false);
    double sigma = 0; for (uint32_t x : sp[0]) sigma += x; sigma /= 100.0; if (sigma < 1) sigma = 1;
    double want = c * (double)g.n / sigma; uint32_t b = 64; while (b < want && b < bmax) b *= 2;
    fprintf(stderr, "auto_b: pilot sigma(top-10)=%.1f n=%u -> b=%u\n", sigma, g.n, b);
    return b;
}

// estimated expected spread of every prefix of `seq` on an (independent) RR sketch: n * covered / #sets
// Need gate: tail-resolution ratio R_tail(k) = b * C_hat_0.05(RIS k-prefix on the validation sketch) / n, k = 1..kmax
static std::vector<double> tail_ratios(const RiskIndex& X, uint32_t b, uint32_t n, uint32_t kmax) {
    std::vector<double> out; const auto& ch = X.chains[X.mean_h()];
    for (uint32_t k = 1; k <= kmax && k <= ch.prefix_sum.size(); ++k) out.push_back((double)b * cvar_from_prefix(ch.prefix_sum[k - 1], 0.05) / (double)n);
    return out;
}
static std::vector<double> prefix_spreads(const GRRSketch& sk, const std::vector<uint32_t>& seq, uint32_t kmax) {
    std::vector<uint8_t> cov(sk.num_sets(), 0); uint64_t c = 0; std::vector<double> out;
    for (uint32_t i = 0; i < kmax && i < seq.size(); ++i) {
        uint32_t u = seq[i];
        for (uint64_t j = sk.inv_off[u]; j < sk.inv_off[u + 1]; ++j) { uint32_t r = sk.inv_rr[j]; if (!cov[r]) { cov[r] = 1; ++c; } }
        out.push_back((double)sk.n * (double)c / (double)sk.num_sets());
    }
    return out;
}

// ------------------------------------------------------------------- ORACLE
// Tier A: exhaustive CVaR optimum on n<=64 graphs.  Selection worlds (A) and
// evaluation worlds (B) are disjoint, so the reported optimum is not in-sample.
static int cmd_oracle(const Args& a) {
    std::string family; Graph g = load_graph(a.get("graph"), family);
    Model model = parse_model(a.get("model", "WC")); Diffusion d(g, model, a.getd("p", 0.1));
    uint32_t k = (uint32_t)a.geti("k", 3);
    auto alphas = dlist(a.get("alphas", "0.05,0.1,0.2,0.5,1.0"));
    uint32_t M = (uint32_t)a.geti("M", 256), b = (uint32_t)a.geti("b", 256), H = (uint32_t)a.geti("H", 8), Kmax = (uint32_t)a.geti("Kmax", 5);
    auto seeds = ilist(a.get("seeds", "0,1,2,3,4"));
    uint32_t RA = (uint32_t)a.geti("RA", 5000), RB = (uint32_t)a.geti("RB", 5000);
    std::string tag = a.get("tag", "oracle");
    FILE* out = fopen(a.get("out", "oracle.csv").c_str(), "a");
    // lite=1: skip the exhaustive search and all baselines; recompute only the index variants (the
    // exhaustive optimum is joined from an earlier full run on the identical evaluation worlds).
    bool lite = a.geti("lite", 0) != 0, hyb = a.geti("hyb", 0) != 0;
    double rho_min = a.getd("rhomin", 0.0), hz = a.getd("hz", 0.0);
    Oracle OA(d, RA, 101), OB(d, RB, 202);
    std::vector<double> optA(alphas.size(), std::numeric_limits<double>::quiet_NaN()); std::vector<std::vector<uint32_t>> Sopt, Sopt1;
    if (!lite) { Sopt = OA.best_sets(k, alphas, optA); std::vector<double> optA1; Sopt1 = OA.best_sets(k, std::vector<double>{1.0}, optA1); }  // CVaR- and mean-optimal sets
    std::vector<double> sp;
    for (long seed : seeds) {
        GRRSketch sk = GRRSketch::build(d, M, b, (uint64_t)seed * 7919 + 1, true);
        GRRSketch ind = GRRSketch::build(d, M, b, (uint64_t)seed * 7919 + 1, false);
        RiskIndex idx = RiskIndex::build(sk, H, Kmax, 0, true);
        RiskIndex rawidx; if (!lite) rawidx = RiskIndex::build(sk, H, Kmax, 0, false);
        std::vector<uint32_t> mean_set = mean_greedy(sk, k);
        std::vector<uint32_t> mean_ind = mean_greedy(ind, k);
        RiskIndex hidx;
        if (hyb) { std::vector<uint32_t> ris_seq = mean_greedy(ind, Kmax); hidx = RiskIndex::build(sk, H, Kmax, 0, true, nullptr, true, hz, &ris_seq); hidx.rho_min = rho_min; hidx.margin_within = true; }
        bool val = a.geti("val", 0) != 0; RiskIndex hvidx;
        RiskIndex hvraw, vidx;
        if (val && hyb) { GRRSketch V = GRRSketch::build(d, M, b, (uint64_t)seed * 7919 + 3, true); hvidx = hidx; hvidx.rematerialize(V);
                          vidx = idx; vidx.rho_min = rho_min; vidx.rematerialize(V);
                          std::vector<uint32_t> ris_seq = mean_greedy(ind, Kmax); hvraw = RiskIndex::build(sk, H, Kmax, 0, false, &hidx.taus, true, 0.0, &ris_seq); hvraw.rematerialize(V); }
        RiskIndex gidx = idx; gidx.rho_min = rho_min;  // the original index with the rho-guard
        // Forward layer (need gate R_tail < theta, cost gate sigma_max) on top of the GRR index (kf = k)
        bool with_fl = a.geti("fl", 0) != 0 && val && hyb; ForwardLayer fl; const char* fl_label = std::isfinite(a.getd("fltheta", 1e300)) ? "RiskIndex-AD" : "RiskIndex-LG"; (void)fl_label;
        if (with_fl) {
            ForwardLayer::Params fp; fp.kf = k; fp.Wsel = (uint32_t)a.geti("flWsel", 1000); fp.Wval = (uint32_t)a.geti("flWval", 1000);
            fp.P = (uint32_t)a.geti("flP", 60); fp.Pscreen = (uint32_t)a.geti("flPscreen", 2000); fp.T = (uint32_t)a.geti("flT", 7);
            fp.alphas = dlist(a.get("flalphas", "0.05,0.1,0.2")); fp.sel_stream = (uint64_t)seed * 7919 + 505; fp.val_stream = (uint64_t)seed * 7919 + 707;
            std::vector<uint32_t> order(g.n); for (uint32_t v = 0; v < g.n; ++v) order[v] = v;
            std::sort(order.begin(), order.end(), [&](uint32_t x, uint32_t y) { return ind.deg(x) != ind.deg(y) ? ind.deg(x) > ind.deg(y) : x < y; });
            std::vector<std::vector<uint32_t>> chs; for (auto& c : hvidx.chains) chs.push_back(c.seeds);
            fp.sigma_max = a.getd("flsmax", 1000);
            fp.theta = a.getd("fltheta", 1e300); auto rt = tail_ratios(hvidx, b, g.n, k); fl = ForwardLayer::build(d, order, chs, (int)hvidx.mean_h(), fp, true, prefix_spreads(ind, mean_greedy(ind, Kmax), k), rt);
            // gate diagnostic (analysis only): tail-resolution ratios of the RIS prefixes and the budgets the layer serves
            fprintf(stderr, "[flgate] %s %s %s seed=%ld kf=%u rtail=", a.get("tag", "").c_str(), a.get("graph").c_str(), model_name(model), seed, fl.kf);
            for (size_t i = 0; i < rt.size(); ++i) { fprintf(stderr, "%s%.3f", i ? "," : "", rt[i]); }
            fprintf(stderr, "\n");
        }
        for (size_t ai = 0; ai < alphas.size(); ++ai) {
            double al = alphas[ai];
            std::map<std::string, std::vector<uint32_t>> sets;
            sets["RIS-IndRR"] = mean_ind;
            sets["RiskIndex"] = idx.query(k, al).seeds;
            sets["Mean-GRR"] = mean_set;
            if (hyb) sets["RiskIndex-H"] = hidx.query(k, al).seeds;
            if (val && hyb) { sets["RiskIndex-HV"] = hvidx.query(k, al).seeds; sets["RiskIndex-HV-raw"] = hvraw.query(k, al).seeds; sets["RiskIndex-V"] = vidx.query(k, al).seeds; }
            if (with_fl) sets[fl_label] = fl.serves(k) ? fl.query(k, al).seeds : hvidx.query(k, al).seeds;
            if (rho_min > 0) sets["RiskIndex-G"] = gidx.query(k, al).seeds;
            if (!lite) {
            sets["Oracle"] = Sopt[ai];
            sets["Oracle-Mean"] = Sopt1[0];
            sets["RiskIndex-raw"] = rawidx.query(k, al).seeds;
            sets["DirectCVaR-GRR"] = direct_cvar_greedy(sk, k, al);
            sets["DirectCVaR-GRR-VM"] = direct_cvar_greedy(sk, k, al, nullptr, true);
            sets["DirectCVaR-IndRR"] = direct_cvar_greedy(ind, k, al);
            sets["SG"] = seq_greedy(sk, k, al, H, false, &mean_set).S;
            }
            for (auto& [name, S] : sets) {
                OB.spreads(S, sp); double meanB = mean_of(sp); double cvB = cvar_unsorted(sp, al);
                OA.spreads(S, sp); double cvA = cvar_unsorted(sp, al);
                fprintf(out, "%s,%s,%s,%u,%lu,%s,%u,%g,%ld,%s,%.6f,%.6f,%.6f,%.6f,%s\n", tag.c_str(), a.get("graph").c_str(), family.c_str(), g.n, (unsigned long)g.m,
                        model_name(model), k, al, seed, name.c_str(), cvB, meanB, cvA, optA[ai], vec_str(S).c_str());
            }
        }
    }
    fclose(out); return 0;
}

// ----------------------------------------------------------------- FIDELITY
// RQ1: estimation fidelity of GRR vs independent RR at equal budget M*b, for
// random seed sets.  Truth from R fresh worlds.
static int cmd_fidelity(const Args& a) {
    std::string family; Graph g = load_graph(a.get("graph"), family);
    Model model = parse_model(a.get("model", "WC")); Diffusion d(g, model, a.getd("p", 0.1));
    uint32_t M = (uint32_t)a.geti("M", 256);
    auto bs = ilist(a.get("bs", "1,4,16,64,256"));
    auto seeds = ilist(a.get("seeds", "0,1,2"));
    uint32_t nsets = (uint32_t)a.geti("nsets", 100), k = (uint32_t)a.geti("k", 10);
    uint32_t R = (uint32_t)a.geti("R", 10000);
    auto alphas = dlist(a.get("alphas", "0.05,0.1,0.2"));
    std::string tag = a.get("tag", "fid"), gname = a.get("name", a.get("graph"));
    FILE* out = fopen(a.get("out", "fidelity.csv").c_str(), "a");
    // seed sets: half uniform random, half degree-biased (so spreads vary across sets)
    Rng rng(4242);
    std::vector<uint32_t> by_deg(g.n); for (uint32_t v = 0; v < g.n; ++v) by_deg[v] = v;
    std::sort(by_deg.begin(), by_deg.end(), [&](uint32_t x, uint32_t y) { return g.outdeg(x) > g.outdeg(y); });
    std::vector<std::vector<uint32_t>> sets; std::vector<std::string> kind;
    for (uint32_t s = 0; s < nsets; ++s) {
        std::set<uint32_t> S;
        uint32_t pool = (s % 2 == 0) ? g.n : std::max<uint32_t>(k * 5, g.n / 50);
        while (S.size() < k) S.insert(by_deg[rng.below(pool)]);
        sets.emplace_back(S.begin(), S.end()); kind.push_back("random");
    }
    // candidate-family sets: prefixes of H=16 threshold chains built on a separate
    // sketch; within one budget they have similar means but different tails,
    // which is exactly the discrimination an online query needs.
    if (a.geti("family", 1)) {
        GRRSketch fsk = GRRSketch::build(d, 256, (uint32_t)a.geti("fb", 2048), 777777, true);
        auto fks = ilist(a.get("fks", "5,10,25,50"));
        RiskIndex fidx = RiskIndex::build(fsk, 16, (uint32_t)fks.back());
        std::set<std::string> seen;
        for (long fk : fks) for (uint32_t h = 0; h < fidx.H; ++h) {
            std::vector<uint32_t> S(fidx.chains[h].seeds.begin(), fidx.chains[h].seeds.begin() + fk);
            std::vector<uint32_t> key = S; std::sort(key.begin(), key.end());
            if (!seen.insert(vec_str(key)).second) continue;
            sets.push_back(S); kind.push_back("family" + std::to_string(fk));
        }
    }
    nsets = (uint32_t)sets.size();
    auto truth = evaluate_sets(d, sets, R, 999, g.m < 20000000);
    std::vector<double> tv(R);
    std::vector<double> true_mean(nsets); std::vector<std::vector<double>> true_cvar(nsets, std::vector<double>(alphas.size()));
    for (uint32_t s = 0; s < nsets; ++s) {
        for (uint32_t w = 0; w < R; ++w) tv[w] = truth[s][w];
        true_mean[s] = mean_of(tv); std::sort(tv.begin(), tv.end());
        for (size_t ai = 0; ai < alphas.size(); ++ai) true_cvar[s][ai] = cvar_sorted(tv, alphas[ai]);
    }
    std::vector<double> X; std::vector<uint8_t> covered; std::vector<uint32_t> touched;
    for (long b : bs) for (int grouped = 1; grouped >= 0; --grouped) for (long seed : seeds) {
        GRRSketch sk = GRRSketch::build(d, M, (uint32_t)b, (uint64_t)seed * 7919 + 11, grouped);
        covered.assign(sk.num_sets(), 0);
        for (uint32_t s = 0; s < nsets; ++s) {
            sk.profile(sets[s], X, covered, touched);
            double em = mean_of(X); std::vector<double> Y = X; double rho = vm_shrink(Y, (double)g.n, (double)b);
            std::sort(X.begin(), X.end()); std::sort(Y.begin(), Y.end());
            std::string line = tag + "," + gname + "," + model_name(model) + "," + std::to_string(M) + "," + std::to_string(b) + "," +
                               (grouped ? "GRR" : "IndRR") + "," + std::to_string(seed) + "," + std::to_string(s) + "," + kind[s] + "," + std::to_string(true_mean[s]) + "," + std::to_string(em);
            for (size_t ai = 0; ai < alphas.size(); ++ai) line += "," + std::to_string(true_cvar[s][ai]) + "," + std::to_string(cvar_sorted(X, alphas[ai])) + "," + std::to_string(cvar_sorted(Y, alphas[ai]));
            fprintf(out, "%s,%.4f,%.4f,%lu\n", line.c_str(), rho, sk.build_sec, (unsigned long)sk.L);
        }
    }
    fclose(out); return 0;
}

// ------------------------------------------------------------------ QUALITY
// Tier B / sensitivity: build sketch + index, answer a (k, alpha) grid with all
// methods, evaluate every selected set on R fresh worlds (common random numbers).
static int cmd_quality(const Args& a) {
    std::string family; Graph g = load_graph(a.get("graph"), family);
    Model model = parse_model(a.get("model", "WC")); Diffusion d(g, model, a.getd("p", 0.01));
    auto Ms = ilist(a.get("Ms", a.get("M", "256"))), bs = ilist(a.get("bs", a.get("b", "2048")));
    auto Hs = ilist(a.get("Hs", a.get("H", "8"))), Ks = ilist(a.get("Kmaxs", a.get("Kmax", "100")));
    double split = a.getd("split", 0.0);  // fraction of worlds used for chain building (0 = in-sample, all worlds)
    auto ks = ilist(a.get("ks", "10,25,50,100")); auto alphas = dlist(a.get("alphas", "0.05,0.1,0.2,0.5,1.0"));
    auto seeds = ilist(a.get("seeds", "0,1,2")); uint32_t R = (uint32_t)a.geti("R", 5000);
    bool with_ind = a.geti("ind", 1) != 0, with_direct = a.geti("direct", 1) != 0, with_nomat = a.geti("nomat", 1) != 0, with_inddirect = a.geti("inddirect", 1) != 0;
    bool vm = a.geti("vm", 1) != 0, with_rawidx = a.geti("rawidx", 0) != 0, with_directvm = a.geti("directvm", 0) != 0;
    bool with_directvm_legacy = a.geti("directvm_legacy", 0) != 0;        // old (invalid) pruning bound, for the audit only
    bool with_sg = a.geti("sg", 0) != 0, with_sgc = a.geti("sgc", 0) != 0;  // Zhou-Tokekar sequential greedy: per-query / cached
    uint32_t Hsg = (uint32_t)a.geti("Hsg", 0);                              // SG grid size (0 = same H as the index)
    std::string tag = a.get("tag", "qual"), gname = a.get("name", a.get("graph"));
    if (split > 0) tag += "-split";
    FILE* out = fopen(a.get("out", "quality.csv").c_str(), "a");
    FILE* sys = fopen(a.get("sysout", "quality_sys.csv").c_str(), "a");
    FILE* spf = a.has("spreads") ? fopen(a.get("spreads").c_str(), "a") : nullptr;
    FILE* cf = fopen(a.get("certout", a.get("out", "quality.csv") + ".cert").c_str(), "a");  // certificates + chosen chain
    bool materialize = g.m < 20000000;
    if (a.has("brule")) bs = {(long)auto_b(d, a.getd("brule", 100), (uint32_t)a.geti("bmax", 16384))};
    double bscale = a.getd("bscale", 1.0);  // dev: rescale the rule's b (e.g. 0.5 with M doubled keeps the budget M*b)
    if (bscale != 1.0) for (long& b : bs) b = std::max<long>(2, (long)std::llround(b * bscale));
    for (long M : Ms) for (long b : bs) for (long H : Hs) for (long Kmax : Ks) for (long seed : seeds) {
        std::string cfg = std::to_string(M) + "," + std::to_string(b) + "," + std::to_string(H) + "," + std::to_string(Kmax) + "," + std::to_string(seed);
        GRRSketch sk = GRRSketch::build(d, (uint32_t)M, (uint32_t)b, (uint64_t)seed * 7919 + 101, true);
        double zmargin = a.getd("z", 0.0); bool store_sq = a.geti("sq", zmargin > 0 ? 1 : 0) != 0;
        RiskIndex idx = RiskIndex::build(sk, (uint32_t)H, (uint32_t)Kmax, split > 0 ? (uint32_t)(split * M) : 0, vm, nullptr, store_sq, zmargin);
        // optional additional index variants (same sketch) for selection-rule ablations: zs=0.5,1,2
        std::vector<std::pair<double, RiskIndex>> zvars;
        for (double zz : dlist(a.get("zs", ""))) zvars.push_back({zz, RiskIndex::build(sk, (uint32_t)H, (uint32_t)Kmax, split > 0 ? (uint32_t)(split * M) : 0, vm, nullptr, true, zz)});
        RiskIndex rawidx; if (with_rawidx) rawidx = RiskIndex::build(sk, (uint32_t)H, (uint32_t)Kmax, 0, false);
        GRRSketch ind; double ind_build = 0;
        if (with_ind) { ind = GRRSketch::build(d, (uint32_t)M, (uint32_t)b, (uint64_t)seed * 7919 + 101, false); ind_build = ind.build_sec; }
        // --- hybrid index: tau = inf chain = RIS greedy on the independent sketch, materialized on GRR ---
        bool with_hyb = a.geti("hyb", 0) != 0 && with_ind;
        RiskIndex hidx; double ris_sec = 0; std::vector<std::pair<double, RiskIndex>> hzvars; std::vector<uint32_t> ris_seq;
        if (with_hyb) {
            Timer tr; ris_seq = mean_greedy(ind, (uint32_t)Kmax); ris_sec = tr.sec();
            hidx = RiskIndex::build(sk, (uint32_t)H, (uint32_t)Kmax, 0, vm, nullptr, true, a.getd("hz", 0.0), &ris_seq);
            hidx.rho_min = a.getd("rhomin", 0.0); hidx.margin_within = a.geti("hzw", 1) != 0;
            for (double zz : dlist(a.get("hzs", ""))) { hzvars.push_back({zz, hidx}); hzvars.back().second.z = zz; }
            // dev ablation: the same margin with the between-world standard error only (encoded as negative z)
            for (double zz : dlist(a.get("hzbs", ""))) { hzvars.push_back({-zz, hidx}); hzvars.back().second.z = zz; hzvars.back().second.margin_within = false; }
            fprintf(stderr, "[hyb] RIS chain %.2fs, hybrid index %.2fs\n", ris_sec, hidx.build_sec);
        }
        // --- validation sketch: an independent GRR sketch of the same shape on which the chains are re-materialized ---
        bool with_val = a.geti("val", 0) != 0;
        RiskIndex vidx, hvidx, hvraw; std::vector<std::pair<double, RiskIndex>> hvzvars; double val_sketch_sec = 0;
        bool with_hvraw = a.geti("hvraw", 0) != 0 && with_hyb;
        if (with_val) {
            GRRSketch V = GRRSketch::build(d, (uint32_t)M, (uint32_t)b, (uint64_t)seed * 7919 + 303, true); val_sketch_sec = V.build_sec;
            vidx = idx; vidx.rematerialize(V);
            if (with_hyb) {
                hvidx = hidx; hvidx.z = a.getd("hvz", 0.0); hvidx.rematerialize(V);
                for (double zz : dlist(a.get("hvzs", ""))) { hvzvars.push_back({zz, hvidx}); hvzvars.back().second.z = zz; }
                // ablation: the same index with raw (not variance-matched) profiles on the validation sketch
                if (with_hvraw) { hvraw = RiskIndex::build(sk, (uint32_t)H, (uint32_t)Kmax, 0, false, &hidx.taus, true, 0.0, &ris_seq); hvraw.rematerialize(V); }
            }
            fprintf(stderr, "[val] validation sketch %.2fs\n", val_sketch_sec);
        }
        // --- Forward layer for under-resolved small-budget tails (RiskIndex-AD; fltheta=inf gives the ungated ablation RiskIndex-LG) ---
        bool with_fl = a.geti("fl", 0) != 0 && with_hyb && with_val; const char* fl_label = std::isfinite(a.getd("fltheta", 1e300)) ? "RiskIndex-AD" : "RiskIndex-LG"; (void)fl_label;
        ForwardLayer fl;
        if (with_fl) {
            ForwardLayer::Params fp; fp.kf = (uint32_t)a.geti("flk", 5); fp.Wsel = (uint32_t)a.geti("flWsel", 1000); fp.Wval = (uint32_t)a.geti("flWval", 1000);
            fp.P = (uint32_t)a.geti("flP", 60); fp.Pscreen = (uint32_t)a.geti("flPscreen", 2000); fp.T = (uint32_t)a.geti("flT", 7);
            fp.alphas = dlist(a.get("flalphas", "0.05,0.1,0.2"));
            fp.sel_stream = (uint64_t)seed * 7919 + 505; fp.val_stream = (uint64_t)seed * 7919 + 707;
            std::vector<uint32_t> order(g.n); for (uint32_t v = 0; v < g.n; ++v) order[v] = v;
            std::sort(order.begin(), order.end(), [&](uint32_t x, uint32_t y) { return ind.deg(x) != ind.deg(y) ? ind.deg(x) > ind.deg(y) : x < y; });
            std::vector<std::vector<uint32_t>> chs; for (auto& c : hvidx.chains) chs.push_back(c.seeds);
            fp.sigma_max = a.getd("flsmax", 1000);
            fp.theta = a.getd("fltheta", 1e300); fl = ForwardLayer::build(d, order, chs, (int)hvidx.mean_h(), fp, materialize, prefix_spreads(ind, ris_seq, fp.kf), tail_ratios(hvidx, (uint32_t)b, g.n, fp.kf));
            fprintf(stderr, "[fl] forward layer kf=%u %.1fs (screen %.1f greedy %.1f val %.1f) pool=%zu cands=%zu\n", fl.kf, fl.build_sec, fl.screen_sec, fl.greedy_sec, fl.val_sec, fl.pool_size, fl.cands.size());
        }
        FILE* hcf = (with_hyb || with_val) ? fopen((a.get("out", "quality.csv") + ".hcert").c_str(), "a") : nullptr;
        // --- collect seed sets for every (k, alpha) -----------------------------
        struct Rec { std::string method; long k; double alpha; double sec; std::vector<uint32_t> S; };
        std::vector<Rec> recs;
        std::vector<double> ri_lat, nomat_lat;
        RiskIndexNoMat nomat(sk, idx);
        uint32_t Hsg_eff = Hsg ? Hsg : (uint32_t)H;
        std::unique_ptr<SeqGreedyCache> sgc; double sgc_build = 0;
        if (with_sgc) { sgc.reset(new SeqGreedyCache(sk, (uint32_t)Kmax, Hsg_eff)); sgc_build = sgc->build_sec; }
        std::unique_ptr<SeqGreedyCache> sgg;  // SG selection rule on RiskIndex's own geometric grid (matched thresholds)
        if (a.geti("sgg", 0)) sgg.reset(new SeqGreedyCache(sk, (uint32_t)Kmax, (uint32_t)H, &idx.taus));
        double rho_min = a.getd("rhomin", 0.0); idx.rho_min = rho_min;
        std::vector<std::pair<double, RiskIndex>> rvars;  // rho-guard ablation at several thresholds (same index, copies are cheap enough)
        for (double rr : dlist(a.get("rhomins", ""))) { rvars.push_back({rr, idx}); rvars.back().second.rho_min = rr; }
        for (long k : ks) {
            if (k > Kmax) continue;
            std::vector<uint32_t> mean_grr(idx.chains[idx.mean_h()].seeds.begin(), idx.chains[idx.mean_h()].seeds.begin() + k);
            std::vector<uint32_t> mean_ind; double mean_ind_sec = 0;
            if (with_ind) { Timer t; mean_ind = mean_greedy(ind, (uint32_t)k); mean_ind_sec = t.sec(); }
            for (double al : alphas) {
                { Timer t; auto ans = idx.query((uint32_t)k, al); double s = t.sec(); ri_lat.push_back(s); recs.push_back({"RiskIndex", k, al, s, ans.seeds});
                  fprintf(cf, "%s,%s,%s,%s,%ld,%g,%u,%.6f,%.4f,%.6f,%.6f,%.6f\n", tag.c_str(), gname.c_str(), model_name(model), cfg.c_str(), k, al, ans.h, ans.score, idx.chains[ans.h].rho[k - 1], cvar_from_prefix(idx.chains[idx.mean_h()].prefix_sum[k - 1], al), ans.se, ans.raw_valid ? ans.raw : std::numeric_limits<double>::quiet_NaN()); }
                for (auto& zv : zvars) { auto ans = zv.second.query((uint32_t)k, al); char nm[64]; snprintf(nm, 64, "RiskIndex-z%g", zv.first); recs.push_back({nm, k, al, 0.0, ans.seeds}); }
                if (with_hyb || with_val) {
                    auto hq = [&](const RiskIndex& X, const std::string& nm) {
                        Timer t; auto ans = X.query((uint32_t)k, al); double s = t.sec(); recs.push_back({nm, k, al, s, ans.seeds});
                        uint32_t hm = X.mean_h();
                        fprintf(hcf, "%s,%s,%s,%s,%s,%ld,%g,%u,%.6f,%.4f,%.6f,%.6f,%.6f,%.6f,%d\n", tag.c_str(), gname.c_str(), model_name(model), cfg.c_str(), nm.c_str(), k, al, ans.h, ans.score,
                                ans.rho, cvar_from_prefix(X.chains[hm].prefix_sum[k - 1], al), ans.se, ans.se_w, ans.raw_valid ? ans.raw : std::numeric_limits<double>::quiet_NaN(), (int)ans.guarded);
                    };
                    if (with_hyb) {
                        hq(hidx, "RiskIndex-H");
                        for (auto& zv : hzvars) { char nm[64]; if (zv.first >= 0) snprintf(nm, 64, "RiskIndex-H-z%g", zv.first); else snprintf(nm, 64, "RiskIndex-H-zb%g", -zv.first); hq(zv.second, nm); }
                    }
                    if (with_val) {
                        hq(vidx, "RiskIndex-V");
                        if (with_hyb) { hq(hvidx, "RiskIndex-HV"); for (auto& zv : hvzvars) { char nm[64]; snprintf(nm, 64, "RiskIndex-HV-z%g", zv.first); hq(zv.second, nm); } }
                        if (with_hvraw) hq(hvraw, "RiskIndex-HV-raw");
                        if (with_fl) {
                            if (fl.serves((uint32_t)k)) {
                                Timer t; auto fa = fl.query((uint32_t)k, al); double s = t.sec(); recs.push_back({fl_label, k, al, s, fa.seeds});
                                double ris_c = cvar_from_prefix(fl.cands[fl.ris_cand].prefix[k - 1], al);
                                fprintf(hcf, "%s,%s,%s,%s,%s,%ld,%g,%d,%.6f,%.4f,%.6f,%.6f,%.6f,%.6f,%d\n", tag.c_str(), gname.c_str(), model_name(model), cfg.c_str(), fl_label, k, al, fa.c, fa.score, 1.0, ris_c, fa.se, 0.0, fa.score, 0);
                            } else hq(hvidx, fl_label);
                        }
                    }
                }
                if (with_nomat) { Timer t; nomat.query((uint32_t)k, al); nomat_lat.push_back(t.sec()); }
                if (with_rawidx) { auto ans = rawidx.query((uint32_t)k, al); recs.push_back({"RiskIndex-raw", k, al, 0.0, ans.seeds}); }
                if (with_direct) { GreedyStats st; auto S = direct_cvar_greedy(sk, (uint32_t)k, al, &st); recs.push_back({"DirectCVaR-GRR", k, al, st.sec, S}); }
                if (with_directvm) { GreedyStats st; auto S = direct_cvar_greedy(sk, (uint32_t)k, al, &st, true); recs.push_back({"DirectCVaR-GRR-VM", k, al, st.sec, S}); }
                if (with_directvm_legacy) { GreedyStats st; auto S = direct_cvar_greedy(sk, (uint32_t)k, al, &st, true, true); recs.push_back({"DirectCVaR-GRR-VM-legacy", k, al, st.sec, S}); }
                if (with_sg) { auto r = seq_greedy(sk, (uint32_t)k, al, Hsg_eff, false, &mean_grr); recs.push_back({"SG", k, al, r.stats.sec, r.S}); }
                if (with_sgc) { Timer t; auto S = sgc->query((uint32_t)k, al); recs.push_back({"SG-cached", k, al, t.sec(), S}); }
                if (sgg) { Timer t; auto S = sgg->query((uint32_t)k, al); recs.push_back({"SG-geogrid", k, al, t.sec(), S}); }
                for (auto& rv : rvars) { auto ans = rv.second.query((uint32_t)k, al); char nm[64]; snprintf(nm, 64, "RiskIndex-rho%g", rv.first); recs.push_back({nm, k, al, 0.0, ans.seeds}); }
                if (with_ind && with_inddirect && (with_direct || a.geti("inddirect", 1) == 2)) { GreedyStats st; auto S = direct_cvar_greedy(ind, (uint32_t)k, al, &st); recs.push_back({"DirectCVaR-IndRR", k, al, st.sec, S}); }
                recs.push_back({"Mean-GRR", k, al, 0.0, mean_grr});
                if (with_ind) recs.push_back({"RIS-IndRR", k, al, mean_ind_sec, mean_ind});
            }
        }
        // --- independent evaluation on fresh worlds ---------------------------
        std::vector<std::vector<uint32_t>> uniq; std::map<std::string, size_t> key2id; std::vector<size_t> rec2id(recs.size());
        for (size_t i = 0; i < recs.size(); ++i) {
            std::vector<uint32_t> S = recs[i].S; std::sort(S.begin(), S.end()); std::string key = vec_str(S);
            auto it = key2id.find(key);
            if (it == key2id.end()) { key2id[key] = uniq.size(); rec2id[i] = uniq.size(); uniq.push_back(recs[i].S); } else rec2id[i] = it->second;
        }
        Timer te; auto sp = evaluate_sets(d, uniq, R, 31337 + seed, materialize); double eval_sec = te.sec();
        std::vector<std::vector<double>> stats(uniq.size());  // [mean, cvar(alphas)...]
        std::vector<double> v(R);
        for (size_t i = 0; i < uniq.size(); ++i) {
            for (uint32_t w = 0; w < R; ++w) v[w] = sp[i][w];
            stats[i].push_back(mean_of(v)); std::sort(v.begin(), v.end());
            for (double al : alphas) stats[i].push_back(cvar_sorted(v, al));
            if (spf) { fprintf(spf, "%s,%s,%s,%s,%zu", tag.c_str(), gname.c_str(), model_name(model), cfg.c_str(), i); for (uint32_t w = 0; w < R; ++w) fprintf(spf, ",%u", sp[i][w]); fprintf(spf, "\n"); }
        }
        for (size_t i = 0; i < recs.size(); ++i) {
            const auto& st = stats[rec2id[i]];
            size_t ai = std::find(alphas.begin(), alphas.end(), recs[i].alpha) - alphas.begin();
            fprintf(out, "%s,%s,%s,%s,%s,%ld,%g,%.6f,%.6f,%.6f,%zu,%s\n", tag.c_str(), gname.c_str(), model_name(model), cfg.c_str(), recs[i].method.c_str(), recs[i].k, recs[i].alpha,
                    recs[i].sec, st[0], st[1 + ai], rec2id[i], vec_str(st).c_str());
        }
        // --- systems record ----------------------------------------------------
        double direct_mean = 0; int nd = 0; for (auto& r : recs) if (r.method == "DirectCVaR-GRR") { direct_mean += r.sec; ++nd; }
        if (nd) direct_mean /= nd;
        fprintf(sys, "%s,%s,%s,%u,%lu,%s,%.6f,%.6f,%.6f,%lu,%lu,%lu,%.3e,%.3e,%.3e,%.3e,%.6f,%.6f,%.1f,%lu,%.3f,%s\n", tag.c_str(), gname.c_str(), model_name(model), g.n, (unsigned long)g.m, cfg.c_str(),
                sk.build_sec, idx.build_sec, ind_build, (unsigned long)sk.L, (unsigned long)sk.bytes(), (unsigned long)idx.bytes(),
                percentile(ri_lat, 0.5), percentile(ri_lat, 0.95), percentile(nomat_lat, 0.5), percentile(nomat_lat, 0.95), direct_mean, eval_sec, peak_rss_mb(),
                (unsigned long)idx.gain_evals, (double)sk.L / sk.num_sets(), vec_str(idx.taus).c_str());
        if (with_sgc) fprintf(sys, "sgc_build,%s,%s,%s,%.6f\n", gname.c_str(), model_name(model), cfg.c_str(), sgc_build);
        if (with_hyb) { double rs = 0; int nr = 0; for (auto& c : hidx.chains) for (double r : c.rho) { rs += r; ++nr; }
            fprintf(sys, "hyb,%s,%s,%s,%.6f,%.6f,%.6f,%lu,%.4f,%s\n", gname.c_str(), model_name(model), cfg.c_str(), ind_build, ris_sec, hidx.build_sec, (unsigned long)hidx.gain_evals, nr ? rs / nr : 0.0, vec_str(hidx.taus).c_str()); }
        if (with_val) fprintf(sys, "val,%s,%s,%s,%.6f,%.6f\n", gname.c_str(), model_name(model), cfg.c_str(), val_sketch_sec, with_hyb ? hvidx.build_sec : vidx.build_sec);
        if (with_fl) fprintf(sys, "fl,%s,%s,%s,%.6f,%.6f,%.6f,%.6f,%zu,%zu,%zu,%u\n", gname.c_str(), model_name(model), cfg.c_str(), fl.build_sec, fl.screen_sec, fl.greedy_sec, fl.val_sec, fl.pool_size, fl.cands.size(), fl.bytes(), fl.kf);
        if (hcf) fclose(hcf);
        { double rs = 0; int nr = 0; for (auto& c : idx.chains) for (double r : c.rho) { rs += r; ++nr; } fprintf(sys, "rho,%s,%s,%s,%.4f\n", gname.c_str(), model_name(model), cfg.c_str(), nr ? rs / nr : 0.0); }
        fflush(out); fflush(sys);
        fprintf(stderr, "[%s %s cfg=%s] sketch %.2fs index %.2fs eval %.2fs direct %.3fs L=%lu\n", gname.c_str(), model_name(model), cfg.c_str(), sk.build_sec, idx.build_sec, eval_sec, direct_mean, (unsigned long)sk.L);
    }
    fclose(out); fclose(sys); fclose(cf); if (spf) fclose(spf); return 0;
}

// -------------------------------------------------------------------- SCALE
// Tier C: offline cost, index size, online latency distribution over a random
// workload, direct per-query cost on the (k, alpha) grid, break-even.
static int cmd_scale(const Args& a) {
    std::string family; Timer tl; Graph g = load_graph(a.get("graph"), family); double load_sec = tl.sec();
    Model model = parse_model(a.get("model", "WC")); Diffusion d(g, model, a.getd("p", 0.01));
    uint32_t M = (uint32_t)a.geti("M", 256), b = (uint32_t)a.geti("b", 2048), H = (uint32_t)a.geti("H", 8), Kmax = (uint32_t)a.geti("Kmax", 100);
    if (a.has("brule")) b = auto_b(d, a.getd("brule", 100), (uint32_t)a.geti("bmax", 16384));
    auto ks = ilist(a.get("ks", "10,25,50,100")); auto alphas = dlist(a.get("alphas", "0.05,0.1,0.2,0.5,1.0"));
    auto seeds = ilist(a.get("seeds", "0")); uint32_t Q = (uint32_t)a.geti("Q", 1000), R = (uint32_t)a.geti("R", 0);
    bool with_direct = a.geti("direct", 1) != 0; uint32_t nomat_q = (uint32_t)a.geti("nomatq", 100);
    std::string tag = a.get("tag", "scale"), gname = a.get("name", a.get("graph"));
    FILE* out = fopen(a.get("out", "scale.csv").c_str(), "a");
    FILE* qf = fopen(a.get("qout", "scale_queries.csv").c_str(), "a");
    bool hyb = a.geti("hyb", 0) != 0; bool store_sq = a.geti("sq", 0) != 0 || hyb;
    double rho_min = a.getd("rhomin", 0.0), hz = a.getd("hz", 0.0);
    FILE* hf = hyb ? fopen(a.get("hout", "scale_hyb.csv").c_str(), "a") : nullptr;
    for (long seed : seeds) {
        double rss0 = peak_rss_mb();
        // hybrid: the tau = inf chain is RIS greedy on an independent RR sketch of the same size, which is
        // built first and freed before the GRR sketch, so the peak memory is that of one sketch.
        std::vector<uint32_t> ris_seq, cov_order; std::vector<double> fl_est; double ind_sec = 0, ris_sec = 0;
        bool with_fl = a.geti("fl", 0) != 0 && hyb && a.geti("val", 0) != 0;
        if (hyb) { GRRSketch ind = GRRSketch::build(d, M, b, (uint64_t)seed * 7919 + 101, false); ind_sec = ind.build_sec;
                   Timer tr; ris_seq = mean_greedy(ind, Kmax); ris_sec = tr.sec();
                   if (with_fl) { fl_est = prefix_spreads(ind, ris_seq, (uint32_t)a.geti("flk", 5)); cov_order.resize(g.n); for (uint32_t v = 0; v < g.n; ++v) cov_order[v] = v;
                                  std::sort(cov_order.begin(), cov_order.end(), [&](uint32_t x, uint32_t y) { return ind.deg(x) != ind.deg(y) ? ind.deg(x) > ind.deg(y) : x < y; }); } }
        GRRSketch sk = GRRSketch::build(d, M, b, (uint64_t)seed * 7919 + 101, true);
        RiskIndex idx = RiskIndex::build(sk, H, Kmax, 0, a.geti("vm", 1) != 0, nullptr, store_sq, hz, hyb ? &ris_seq : nullptr);
        idx.rho_min = rho_min; idx.margin_within = true;
        if (hyb) idx.build_sec += ind_sec + ris_sec;  // the independent sketch and the RIS chain are charged to the index
        double val_sec = 0;
        if (a.geti("val", 0) != 0) {  // validation sketch (freed after re-materialization); charged to the index
            GRRSketch V = GRRSketch::build(d, M, b, (uint64_t)seed * 7919 + 303, true); val_sec = V.build_sec;
            idx.rematerialize(V); idx.build_sec += val_sec;
        }
        // Forward layer (charged to the index); small-budget lookup latency
        ForwardLayer fl; double fl_sec = 0; std::vector<double> fllat;
        if (with_fl) {
            ForwardLayer::Params fp; fp.kf = (uint32_t)a.geti("flk", 5); fp.Wsel = (uint32_t)a.geti("flWsel", 1000); fp.Wval = (uint32_t)a.geti("flWval", 1000);
            fp.P = (uint32_t)a.geti("flP", 60); fp.Pscreen = (uint32_t)a.geti("flPscreen", 2000); fp.T = (uint32_t)a.geti("flT", 7);
            fp.alphas = dlist(a.get("flalphas", "0.05,0.1,0.2")); fp.sel_stream = (uint64_t)seed * 7919 + 505; fp.val_stream = (uint64_t)seed * 7919 + 707;
            std::vector<std::vector<uint32_t>> chs; for (auto& c : idx.chains) chs.push_back(c.seeds);
            fp.sigma_max = a.getd("flsmax", 1000);
            fp.theta = a.getd("fltheta", 1e300); fl = ForwardLayer::build(d, cov_order, chs, (int)idx.mean_h(), fp, g.m < 20000000, fl_est, tail_ratios(idx, (uint32_t)b, g.n, fp.kf)); fl_sec = fl.build_sec; idx.build_sec += fl_sec;
            Rng fr(777 + seed);
            for (uint32_t q = 0; q < Q && fl.kf > 0; ++q) {
                uint32_t k = 1 + fr.below(fl.kf); double al = alphas[fr.below((uint32_t)alphas.size())];
                auto t0 = std::chrono::steady_clock::now(); auto ans = fl.query(k, al); auto t1 = std::chrono::steady_clock::now();
                fllat.push_back(std::chrono::duration<double>(t1 - t0).count()); asm volatile("" : : "r"(ans.seeds.data()) : "memory");
            }
            fprintf(stderr, "[fl] forward layer kf=%u %.1fs (screen %.1f greedy %.1f val %.1f) pool=%zu small-k lookup P50 %.2e\n", fl.kf, fl.build_sec, fl.screen_sec, fl.greedy_sec, fl.val_sec, fl.pool_size, fllat.empty() ? 0.0 : percentile(fllat, 0.5));
        }
        double rss1 = peak_rss_mb();
        // online workload: Q random (k, alpha) from the grid
        Rng wr(555 + seed); std::vector<double> lat; lat.reserve(Q);
        Timer tw;
        for (uint32_t q = 0; q < Q; ++q) {
            uint32_t k = (uint32_t)ks[wr.below((uint32_t)ks.size())]; double al = alphas[wr.below((uint32_t)alphas.size())];
            auto t0 = std::chrono::steady_clock::now(); auto ans = idx.query(k, al); auto t1 = std::chrono::steady_clock::now();
            lat.push_back(std::chrono::duration<double>(t1 - t0).count());
            asm volatile("" : : "r"(ans.seeds.data()) : "memory");
        }
        double workload_sec = tw.sec();
        // no-materialization lookup latency
        RiskIndexNoMat nomat(sk, idx); std::vector<double> nlat;
        for (uint32_t q = 0; q < nomat_q; ++q) {
            uint32_t k = (uint32_t)ks[wr.below((uint32_t)ks.size())]; double al = alphas[wr.below((uint32_t)alphas.size())];
            Timer t; nomat.query(k, al); nlat.push_back(t.sec());
        }
        // direct per-query cost on the grid (+ optional quality on R fresh worlds)
        std::vector<double> dlat; std::vector<std::vector<uint32_t>> sets; std::vector<std::string> labels;
        for (long k : ks) for (double al : alphas) {
            auto ans = idx.query((uint32_t)k, al);
            sets.push_back(ans.seeds); labels.push_back("RiskIndex," + std::to_string(k) + "," + std::to_string(al));
            if (hyb) { sets.push_back(std::vector<uint32_t>(ris_seq.begin(), ris_seq.begin() + k)); labels.push_back("RIS-IndRR," + std::to_string(k) + "," + std::to_string(al));
                       fprintf(hf, "%s,%s,%s,%ld,%ld,%g,%u,%.6f,%.4f,%.6f,%.6f,%d\n", tag.c_str(), gname.c_str(), model_name(model), seed, k, al, ans.h, ans.score, ans.rho, ans.se, ans.se_w, (int)ans.guarded); }
            if (with_direct) { GreedyStats st; auto S = direct_cvar_greedy(sk, (uint32_t)k, al, &st); dlat.push_back(st.sec); sets.push_back(S); labels.push_back("DirectCVaR-GRR," + std::to_string(k) + "," + std::to_string(al));
                fprintf(qf, "%s,%s,%s,%ld,%g,%ld,%.6f,%lu,%lu\n", tag.c_str(), gname.c_str(), model_name(model), k, al, seed, st.sec, (unsigned long)st.cvar_evals, (unsigned long)st.mean_evals); }
        }
        double direct_mean = dlat.empty() ? 0 : mean_of(dlat);
        double lookup_mean = mean_of(lat);
        double qstar = direct_mean > lookup_mean ? (sk.build_sec + idx.build_sec) / (direct_mean - lookup_mean) : -1;
        fprintf(out, "%s,%s,%s,%u,%lu,%u,%u,%u,%u,%ld,%.3f,%.4f,%.4f,%lu,%lu,%lu,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e,%.4f,%.3e,%.3e,%.3e,%.3f,%.1f,%.1f,%.3f\n", tag.c_str(), gname.c_str(), model_name(model), g.n, (unsigned long)g.m, M, b, H, Kmax, seed,
                load_sec, sk.build_sec, idx.build_sec, (unsigned long)sk.L, (unsigned long)sk.bytes(), (unsigned long)idx.bytes(),
                percentile(lat, 0.5), percentile(lat, 0.95), lookup_mean, workload_sec, percentile(nlat, 0.5), percentile(nlat, 0.95),
                direct_mean, percentile(dlat, 0.5), percentile(dlat, 0.95), dlat.empty() ? 0 : *std::max_element(dlat.begin(), dlat.end()), qstar, rss0, rss1, (double)sk.L / sk.num_sets());
        if (hyb) { fprintf(hf, "build,%s,%s,%ld,%.4f,%.4f,%.4f,%.4f,%lu,%.4f,%.1f\n", gname.c_str(), model_name(model), seed, ind_sec, ris_sec, sk.build_sec, idx.build_sec, (unsigned long)idx.bytes(), val_sec, rss1); fflush(hf); }
        if (with_fl) { fprintf(hf, "fl,%s,%s,%ld,%.4f,%.4f,%.4f,%.4f,%zu,%zu,%.3e,%.3e,%u\n", gname.c_str(), model_name(model), seed, fl.build_sec, fl.screen_sec, fl.greedy_sec, fl.val_sec, fl.pool_size, fl.bytes(), fllat.empty() ? 0.0 : percentile(fllat, 0.5), fllat.empty() ? 0.0 : percentile(fllat, 0.95), fl.kf); fflush(hf); }
        fflush(out); fflush(qf);
        fprintf(stderr, "[%s] n=%u m=%lu sketch %.2fs index %.2fs lookup P50 %.2e P95 %.2e direct mean %.3fs Q*=%.2f RSS %.0fMB\n", gname.c_str(), g.n, (unsigned long)g.m, sk.build_sec, idx.build_sec, percentile(lat, 0.5), percentile(lat, 0.95), direct_mean, qstar, rss1);
        if (R > 0) {
            auto sp = evaluate_sets(d, sets, R, 31337 + seed, g.m < 20000000);
            FILE* ef = fopen(a.get("eout", "scale_quality.csv").c_str(), "a");
            std::vector<double> v(R);
            for (size_t i = 0; i < sets.size(); ++i) {
                for (uint32_t w = 0; w < R; ++w) v[w] = sp[i][w];
                double mean = mean_of(v); std::sort(v.begin(), v.end());
                std::string s; for (double al : alphas) s += "," + std::to_string(cvar_sorted(v, al));
                fprintf(ef, "%s,%s,%s,%ld,%s,%.4f%s\n", tag.c_str(), gname.c_str(), model_name(model), seed, labels[i].c_str(), mean, s.c_str());
            }
            fclose(ef);
        }
    }
    fclose(out); fclose(qf); if (hf) fclose(hf); return 0;
}

// ------------------------------------------------------------------- GRAPH INFO
// spectral radius of the adjacency matrix by power iteration (mean-field IC
// epidemic threshold p_c = 1/lambda_1)
static double spectral_radius(const Graph& g, int iters = 200) {
    std::vector<double> x(g.n, 1.0), y(g.n);
    double lam = 0;
    for (int it = 0; it < iters; ++it) {
        double norm = 0;
        for (uint32_t v = 0; v < g.n; ++v) { double s = 0; for (uint64_t i = g.rev_off[v]; i < g.rev_off[v + 1]; ++i) s += x[g.rev_src[i]]; y[v] = s; norm += s * s; }
        norm = std::sqrt(norm); if (norm == 0) return 0;
        lam = norm; for (uint32_t v = 0; v < g.n; ++v) x[v] = y[v] / norm;
    }
    return lam;
}
static int cmd_info(const Args& a) {
    std::string family; Graph g = load_graph(a.get("graph"), family);
    uint64_t maxd = 0, maxin = 0; uint32_t iso = 0; for (uint32_t v = 0; v < g.n; ++v) { maxd = std::max<uint64_t>(maxd, g.outdeg(v)); maxin = std::max<uint64_t>(maxin, g.indeg(v)); if (g.outdeg(v) + g.indeg(v) == 0) ++iso; }
    double lam = a.geti("spectral", 1) ? spectral_radius(g) : 0;
    printf("%s,%u,%lu,%.3f,%lu,%lu,%u,%.3f,%.5f\n", a.get("name", a.get("graph")).c_str(), g.n, (unsigned long)g.m, (double)g.m / g.n, (unsigned long)maxd, (unsigned long)maxin, iso, lam, lam > 0 ? 1.0 / lam : 0);
    return 0;
}

// ---------------------------------------------------------------- DIST
// Spread distribution of given seed sets on R fresh worlds (for figures).
static int cmd_dist(const Args& a) {
    std::string family; Graph g = load_graph(a.get("graph"), family);
    Model model = parse_model(a.get("model", "WC")); Diffusion d(g, model, a.getd("p", 0.1));
    std::vector<std::vector<uint32_t>> sets;
    for (auto& part : split(a.get("sets"), ';')) { std::vector<uint32_t> S; for (auto& t : split(part, ' ')) if (!t.empty()) S.push_back((uint32_t)std::stol(t)); sets.push_back(S); }
    uint32_t R = (uint32_t)a.geti("R", 10000);
    auto sp = evaluate_sets(d, sets, R, a.geti("seed", 909), g.m < 20000000);
    for (size_t i = 0; i < sets.size(); ++i) { printf("%s", vec_str(sets[i]).c_str()); for (uint32_t w = 0; w < R; ++w) printf(",%u", sp[i][w]); printf("\n"); }
    return 0;
}


// ---------------------------------------------------------------- PAIR
// Exhaustive search on a tiny graph (n <= 64) for two k-sets with (nearly) equal
// expected spread and maximally different lower-tail CVaR_alpha (Figure 1a).
static int cmd_pair(const Args& a) {
    std::string family; Graph g = load_graph(a.get("graph"), family);
    Model model = parse_model(a.get("model", "IC")); Diffusion d(g, model, a.getd("p", 0.2));
    uint32_t k = (uint32_t)a.geti("k", 3), R = (uint32_t)a.geti("R", 4000); double alpha = a.getd("alpha", 0.1), tol = a.getd("tol", 0.01);
    Oracle O(d, R, a.geti("seed", 7));
    std::vector<std::vector<uint32_t>> sets; std::vector<double> mean, cv;
    std::vector<uint32_t> S(k); for (uint32_t i = 0; i < k; ++i) S[i] = i;
    std::vector<double> sp;
    while (true) {
        O.spreads(S, sp); std::sort(sp.begin(), sp.end());
        sets.push_back(S); mean.push_back(mean_of(sp)); cv.push_back(cvar_sorted(sp, alpha));
        int i = (int)k - 1; while (i >= 0 && S[i] == g.n - k + i) --i; if (i < 0) break;
        S[i]++; for (uint32_t j = i + 1; j < k; ++j) S[j] = S[j - 1] + 1;
    }
    // only consider sets in the top quartile of expected spread (realistic candidates)
    std::vector<double> ms = mean; std::sort(ms.begin(), ms.end()); double floor_ = ms[(size_t)(0.75 * ms.size())];
    std::vector<size_t> idx; for (size_t i = 0; i < sets.size(); ++i) if (mean[i] >= floor_) idx.push_back(i);
    std::sort(idx.begin(), idx.end(), [&](size_t x, size_t y) { return mean[x] < mean[y]; });
    double best = 0; size_t bi = 0, bj = 0;
    for (size_t x = 0; x < idx.size(); ++x)
        for (size_t y = x + 1; y < idx.size() && mean[idx[y]] - mean[idx[x]] <= tol * mean[idx[x]]; ++y) {
            double r = std::max(cv[idx[x]], cv[idx[y]]) / std::max(1e-9, std::min(cv[idx[x]], cv[idx[y]]));
            if (r > best) { best = r; bi = idx[x]; bj = idx[y]; }
        }
    fprintf(stderr, "pair: %s (mean %.2f cvar %.2f) vs %s (mean %.2f cvar %.2f) ratio %.3f\n", vec_str(sets[bi]).c_str(), mean[bi], cv[bi], vec_str(sets[bj]).c_str(), mean[bj], cv[bj], best);
    printf("%s;%s\n", vec_str(sets[bi]).c_str(), vec_str(sets[bj]).c_str());
    return 0;
}

// Self-test: (1) lazy Direct-CVaR greedy (raw and variance-matched objective, valid bounds) must select the same
// sequence as exhaustive greedy without pruning; the legacy 1/alpha bound under VM is reported, not asserted.
// (2) cvar_from_prefix == cvar_unsorted on stored profiles; (3) raw certificate recovered from (rho, mean, C_vm)
// equals the raw empirical CVaR of the uncorrected profile.
static int cmd_selftest(Args& a) {
    int fails = 0, legacy_diff = 0, total = 0;
    for (int inst = 1; inst <= (int)a.geti("inst", 6); ++inst) {
        for (std::string spec : {"er:40:4:" + std::to_string(inst), "ba:48:2:" + std::to_string(inst), "sbm:40:4:0.5:0.03:" + std::to_string(inst)}) {
            std::string fam; Graph g = load_graph(spec, fam);
            for (Model model : {Model::WC, Model::IC, Model::LT}) {
                Diffusion d(g, model, model == Model::IC ? 0.2 : 0.0);
                GRRSketch sk = GRRSketch::build(d, (uint32_t)a.geti("M", 96), (uint32_t)a.geti("b", 24), 1000 + inst, true);
                for (double al : {0.05, 0.1, 0.25, 0.5, 1.0}) {
                    uint32_t k = 4;
                    for (bool vm : {false, true}) {
                        auto lazy = direct_cvar_greedy(sk, k, al, nullptr, vm, false, false);
                        auto full = direct_cvar_greedy(sk, k, al, nullptr, vm, false, true);
                        ++total;
                        auto obj = [&](const std::vector<uint32_t>& S) {  // objective value of a set on the sketch (ties may change the sequence, not the value)
                            std::vector<double> X; std::vector<uint8_t> cov(sk.num_sets(), 0); std::vector<uint32_t> tch;
                            sk.profile(S, X, cov, tch); if (vm) vm_shrink(X, (double)sk.n, (double)sk.b); return cvar_unsorted(X, al); };
                        if (lazy != full && std::fabs(obj(lazy) - obj(full)) > 1e-9 * (1 + std::fabs(obj(full)))) { ++fails; fprintf(stderr, "MISMATCH %s %s alpha=%g vm=%d\n", spec.c_str(), model_name(model), al, (int)vm); }
                        if (vm) { auto leg = direct_cvar_greedy(sk, k, al, nullptr, true, true, false); if (leg != full) ++legacy_diff; }
                    }
                }
                // index consistency
                RiskIndex idx = RiskIndex::build(sk, 4, 6, 0, true, nullptr, true, 0.0);
                RiskIndex raw = RiskIndex::build(sk, 4, 6, 0, false, &idx.taus, true, 0.0);
                std::vector<double> X; std::vector<uint8_t> cov(sk.num_sets(), 0); std::vector<uint32_t> touched;
                for (uint32_t h = 0; h < 4; ++h) for (uint32_t k = 1; k <= 6; ++k) for (double al : {0.05, 0.2, 1.0}) {
                    std::vector<uint32_t> S(idx.chains[h].seeds.begin(), idx.chains[h].seeds.begin() + k);
                    sk.profile(S, X, cov, touched);
                    double craw = cvar_unsorted(X, al);
                    double craw_idx = cvar_from_prefix(raw.chains[h].prefix_sum[k - 1], al);
                    double cvm = cvar_from_prefix(idx.chains[h].prefix_sum[k - 1], al);
                    bool valid; double rec = idx.raw_cert(h, k, cvm, &valid);
                    if (std::fabs(craw - craw_idx) > 1e-6 * (1 + std::fabs(craw))) { ++fails; fprintf(stderr, "PREFIX MISMATCH\n"); }
                    if (valid && std::fabs(rec - craw) > 1e-6 * (1 + std::fabs(craw))) { ++fails; fprintf(stderr, "RAW RECOVERY MISMATCH %g vs %g rho=%g\n", rec, craw, idx.chains[h].rho[k - 1]); }
                    if (!valid && idx.chains[h].rho[k - 1] > 0) { ++fails; fprintf(stderr, "RAW VALIDITY MISMATCH\n"); }
                    // within-world standard error from (P, Q, mu, rho, n, b) vs. a direct computation on the raw profile
                    if (idx.chains[h].rho[k - 1] > 0) {
                        std::vector<double> Y = X; double mu = mean_of(Y); double rho = vm_shrink(Y, (double)sk.n, (double)sk.b);
                        std::vector<size_t> ord(Y.size()); for (size_t i = 0; i < ord.size(); ++i) ord[i] = i;
                        std::stable_sort(ord.begin(), ord.end(), [&](size_t x, size_t y) { return Y[x] < Y[y]; });
                        size_t r = (size_t)std::ceil(al * Y.size() - 1e-12); if (r == 0) r = 1;
                        double sz = 0, szz = 0, sv = 0;
                        for (size_t i = 0; i < r; ++i) { double z = Y[ord[i]], x = X[ord[i]]; sz += z; szz += z * z; sv += x * (sk.n - x) / (sk.b - 1.0); }
                        double m1 = sz / r, var = szz / r - m1 * m1; if (var < 0) var = 0;
                        double direct = std::sqrt(var / r + rho * rho * sv / ((double)r * r));
                        double viaQ = tail_se_within(idx.chains[h].prefix_sum[k - 1], idx.chains[h].prefix_sq[k - 1], al, mu, rho, (double)sk.n, (double)sk.b);
                        if (std::fabs(direct - viaQ) > 1e-6 * (1 + direct)) { ++fails; fprintf(stderr, "SE_W MISMATCH %g vs %g\n", direct, viaQ); }
                    }
                }
                // hybrid: materializing the mean chain's own seed sequence reproduces its stored entries exactly
                {
                    uint32_t hm = idx.mean_h();
                    RiskIndex hy = RiskIndex::build(sk, 4, 6, 0, true, nullptr, true, 0.0, &idx.chains[hm].seeds);
                    for (uint32_t k = 1; k <= 6; ++k) {
                        if (hy.chains[hy.mean_h()].prefix_sum[k - 1] != idx.chains[hm].prefix_sum[k - 1] || hy.chains[hy.mean_h()].rho[k - 1] != idx.chains[hm].rho[k - 1]) { ++fails; fprintf(stderr, "HYBRID MATERIALIZATION MISMATCH\n"); }
                        for (double al : {0.05, 0.2, 1.0}) { auto x1 = hy.query(k, al), x2 = idx.query(k, al); if (x1.seeds != x2.seeds || x1.score != x2.score) { ++fails; fprintf(stderr, "HYBRID QUERY MISMATCH\n"); } }
                    }
                    ++total;
                    // validation: re-materializing on the same sketch reproduces every stored entry
                    RiskIndex rv = idx; rv.rematerialize(sk);
                    for (uint32_t h = 0; h < 4; ++h) for (uint32_t k = 1; k <= 6; ++k)
                        if (rv.chains[h].prefix_sum[k - 1] != idx.chains[h].prefix_sum[k - 1] || rv.chains[h].prefix_sq[k - 1] != idx.chains[h].prefix_sq[k - 1] || rv.chains[h].rho[k - 1] != idx.chains[h].rho[k - 1]) { ++fails; fprintf(stderr, "REMATERIALIZATION MISMATCH\n"); }
                    ++total;
                }
            }
        }
    }
    // query() / query_score() / budget_query() agree: same chain and certificate for every (k, alpha), guard included,
    // and the budget query returns the smallest k whose certificate reaches the target.
    {
        int qf = 0, qt = 0;
        for (int inst = 0; inst < 3; ++inst) {
            Graph g = Graph::erdos_renyi(300, 6.0, 900 + inst); Diffusion d(g, Model::IC, 0.08);
            GRRSketch sk = GRRSketch::build(d, 64, 32, 4242 + inst, true);
            RiskIndex ix = RiskIndex::build(sk, 4, 8, 0, true, nullptr, true, 0.0, nullptr); ix.rho_min = 0.5;
            for (uint32_t k = 1; k <= 8; ++k) for (double al : {0.05, 0.1, 0.2, 0.5, 1.0}) {
                auto q = ix.query(k, al); auto qs = ix.query_score(k, al); ++qt;
                if (q.h != qs.first || q.score != qs.second) { ++qf; fprintf(stderr, "QUERY/QUERY_SCORE MISMATCH k=%u a=%g\n", k, al); }
                auto bq = ix.budget_query(al, q.score); ++qt;
                bool ok = bq.first >= 1 && bq.first <= k && ix.query(bq.first, al).score >= q.score && ix.query(bq.first, al).h == bq.second;
                for (uint32_t kk = 1; ok && kk < bq.first; ++kk) if (ix.query(kk, al).score >= q.score) ok = false;
                if (!ok) { ++qf; fprintf(stderr, "BUDGET_QUERY MISMATCH k=%u a=%g -> %u\n", k, al, bq.first); }
            }
        }
        total += qt; fails += qf;
        printf("selftest: %d query/query_score/budget_query checks, %d failures\n", qt, qf);
    }
    // Routed interfaces with a forward layer: for forward-served (k <= kf) and GRR-served (k > kf) budgets,
    // adaptive_query_score equals the certificate of adaptive_query, and adaptive_budget_query returns the smallest k
    // whose routed certificate reaches the target.
    {
        int qf = 0, qt = 0, nfwd = 0;
        for (int inst = 0; inst < 2; ++inst) {
            Graph g = Graph::erdos_renyi(400, 5.0, 1200 + inst); Diffusion d(g, Model::IC, 0.08);
            GRRSketch sk = GRRSketch::build(d, 64, 32, 5252 + inst, true), ind = GRRSketch::build(d, 1, 64 * 32, 5353 + inst, false);
            std::vector<uint32_t> ris = mean_greedy(ind, 8);
            RiskIndex ix = RiskIndex::build(sk, 4, 8, 0, true, nullptr, true, 0.0, &ris); ix.rho_min = 0.5;
            GRRSketch V = GRRSketch::build(d, 64, 32, 5454 + inst, true); ix.rematerialize(V);
            ForwardLayer::Params fp; fp.kf = 3; fp.Wsel = 200; fp.Wval = 200; fp.P = 10; fp.Pscreen = 100; fp.T = 3; fp.theta = 1e300;
            std::vector<uint32_t> order(g.n); for (uint32_t v = 0; v < g.n; ++v) order[v] = v;
            std::sort(order.begin(), order.end(), [&](uint32_t x, uint32_t y) { return ind.deg(x) != ind.deg(y) ? ind.deg(x) > ind.deg(y) : x < y; });
            std::vector<std::vector<uint32_t>> chs; for (auto& c : ix.chains) chs.push_back(c.seeds);
            std::vector<double> es(8, 0.0);  // cost gate open
            ForwardLayer fl = ForwardLayer::build(d, order, chs, (int)ix.mean_h(), fp, true, es);
            if (fl.kf != 3) { ++qf; fprintf(stderr, "FORWARD LAYER kf=%u (expected 3)\n", fl.kf); }
            for (uint32_t k = 1; k <= 8; ++k) for (double al : {0.05, 0.1, 0.2, 0.5, 1.0}) {
                auto q = adaptive_query(ix, fl, k, al); double qs = adaptive_query_score(ix, fl, k, al); ++qt; if (q.forward) ++nfwd;
                if (q.score != qs || q.forward != fl.serves(k) || (q.forward && (q.rho != 1.0 || q.guarded))) { ++qf; fprintf(stderr, "ADAPTIVE QUERY MISMATCH k=%u a=%g\n", k, al); }
                uint32_t bk = adaptive_budget_query(ix, fl, al, q.score); ++qt;
                bool ok = bk >= 1 && bk <= k && adaptive_query_score(ix, fl, bk, al) >= q.score;
                for (uint32_t kk = 1; ok && kk < bk; ++kk) if (adaptive_query_score(ix, fl, kk, al) >= q.score) ok = false;
                if (!ok) { ++qf; fprintf(stderr, "ADAPTIVE BUDGET_QUERY MISMATCH k=%u a=%g -> %u\n", k, al, bk); }
            }
        }
        total += qt; fails += qf;
        printf("selftest: %d routed query/query_score/budget_query checks (%d forward-served), %d failures\n", qt, nfwd, qf);
    }
    printf("selftest: %d greedy comparisons, %d failures; legacy 1/alpha bound under VM differs from exhaustive in %d cases\n", total, fails, legacy_diff);
    return fails ? 1 : 0;
}

// ---------------------------------------------------------------- HEADROOM
// Headroom diagnosis (development seeds only).  How much tail improvement over RIS exists at small k, and how
// much of it does RiskIndex capture?  Oracle-CVaR: greedy on the exact empirical CVaR of forward-simulated
// spreads on Wsel fresh worlds over a candidate pool (top-P nodes by RR coverage); Oracle-mean: the same greedy
// on the mean (control).  All sets are compared on Weval further fresh worlds (a disjoint stream) together with
// RIS (ind. RR, Mb sets) and RiskIndex-HV built exactly as in the Tier B runs.
static int cmd_headroom(const Args& a) {
    std::string family; Graph g = load_graph(a.get("graph"), family);
    Model model = parse_model(a.get("model", "WC")); Diffusion d(g, model, a.getd("p", 0.01));
    uint32_t P = (uint32_t)a.geti("P", 100), Wsel = (uint32_t)a.geti("Wsel", 1000), Weval = (uint32_t)a.geti("Weval", 3000);
    auto ks = ilist(a.get("ks", "1,2,3,5")); auto alphas = dlist(a.get("alphas", "0.05,0.1,0.2"));
    long seed = a.geti("seed", 100); uint32_t M = (uint32_t)a.geti("M", 1024), H = (uint32_t)a.geti("H", 8), Kmax = (uint32_t)a.geti("Kmax", 100);
    uint32_t b = a.has("brule") ? auto_b(d, a.getd("brule", 100), (uint32_t)a.geti("bmax", 65536)) : (uint32_t)a.geti("b", 1024);
    std::string gname = a.get("name", a.get("graph"));
    FILE* out = fopen(a.get("out", "headroom.csv").c_str(), "a");
    bool materialize = g.m < 20000000;
    uint32_t kmax = 0; for (long k : ks) kmax = std::max<uint32_t>(kmax, (uint32_t)k);
    Timer T;
    GRRSketch sk = GRRSketch::build(d, M, b, (uint64_t)seed * 7919 + 101, true);
    GRRSketch ind = GRRSketch::build(d, M, b, (uint64_t)seed * 7919 + 101, false);
    std::vector<uint32_t> ris_seq = mean_greedy(ind, Kmax);
    RiskIndex hidx = RiskIndex::build(sk, H, Kmax, 0, true, nullptr, true, 0.0, &ris_seq);
    hidx.rho_min = 0.5; hidx.margin_within = true;
    GRRSketch V = GRRSketch::build(d, M, b, (uint64_t)seed * 7919 + 303, true);
    RiskIndex hv = hidx; hv.z = 0.0; hv.rematerialize(V);
    double t_index = T.sec();
    // Forward layer exactly as in the Tier B runs (same streams for the same seed)
    ForwardLayer fl; bool with_fl = a.geti("fl", 1) != 0; const char* fl_label = std::isfinite(a.getd("fltheta", 1e300)) ? "RiskIndex-AD" : "RiskIndex-LG"; (void)fl_label;
    if (with_fl) {
        ForwardLayer::Params fp; fp.kf = (uint32_t)a.geti("flk", 5); fp.sigma_max = a.getd("flsmax", 1000);
        fp.sel_stream = (uint64_t)seed * 7919 + 505; fp.val_stream = (uint64_t)seed * 7919 + 707;
        std::vector<uint32_t> ord(g.n); for (uint32_t v = 0; v < g.n; ++v) ord[v] = v;
        std::sort(ord.begin(), ord.end(), [&](uint32_t x, uint32_t y) { return ind.deg(x) != ind.deg(y) ? ind.deg(x) > ind.deg(y) : x < y; });
        std::vector<std::vector<uint32_t>> chs; for (auto& c : hv.chains) chs.push_back(c.seeds);
        fp.theta = a.getd("fltheta", 1e300); fl = ForwardLayer::build(d, ord, chs, (int)hv.mean_h(), fp, materialize, prefix_spreads(ind, ris_seq, fp.kf), tail_ratios(hv, (uint32_t)b, g.n, fp.kf));
    }
    // candidate pool (screen-then-verify): exact single-seed mean and CVaR of the Pscreen best-covered nodes on the
    // selection worlds; pool = top-P by exact mean, top-P by exact CVaR_alpha for every alpha, the RIS prefix and
    // every seed that RiskIndex-HV returns for k <= kmax.
    uint32_t Pscreen = (uint32_t)a.geti("Pscreen", 2000);
    std::vector<uint32_t> order(g.n); for (uint32_t v = 0; v < g.n; ++v) order[v] = v;
    std::sort(order.begin(), order.end(), [&](uint32_t x, uint32_t y) { return ind.deg(x) > ind.deg(y); });
    uint32_t NS = std::min<uint32_t>(Pscreen, g.n);
    std::vector<std::vector<uint32_t>> singles; for (uint32_t i = 0; i < NS; ++i) singles.push_back({order[i]});
    auto ssp = evaluate_sets(d, singles, Wsel, 555000 + (uint64_t)seed, materialize);
    std::vector<std::vector<double>> sstat(NS);  // [mean, cvar(alphas)...]
    { std::vector<double> v(Wsel);
      for (uint32_t i = 0; i < NS; ++i) { for (uint32_t w = 0; w < Wsel; ++w) v[w] = ssp[i][w]; sstat[i].push_back(mean_of(v)); std::sort(v.begin(), v.end()); for (double al : alphas) sstat[i].push_back(cvar_sorted(v, al)); } }
    std::set<uint32_t> poolset;
    for (size_t c = 0; c <= alphas.size(); ++c) {
        std::vector<uint32_t> idx(NS); for (uint32_t i = 0; i < NS; ++i) idx[i] = i;
        std::sort(idx.begin(), idx.end(), [&](uint32_t x, uint32_t y) { return sstat[x][c] > sstat[y][c]; });
        for (uint32_t i = 0; i < std::min<uint32_t>(P, NS); ++i) poolset.insert(order[idx[i]]);
    }
    for (uint32_t i = 0; i < kmax && i < ris_seq.size(); ++i) poolset.insert(ris_seq[i]);
    for (long k : ks) for (double al : alphas) for (uint32_t u : hv.query((uint32_t)k, al).seeds) poolset.insert(u);
    std::vector<uint32_t> pool(poolset.begin(), poolset.end());
    // forward-simulated references on Wsel worlds (stream 555000+seed), over the pool:
    //   objective 0 = mean, 1..A = CVaR_alpha (myopic CVaR greedy), A+1..A+T = truncated spread F_tau (exact basis),
    // tau grid: T values geometric between the 5% quantile of the RIS 1-prefix and the 20% quantile of the RIS
    // kmax-prefix (exact spreads on the selection worlds).
    uint32_t NT = (uint32_t)a.geti("T", 7);
    std::vector<std::vector<uint32_t>> risk_sets; for (uint32_t k = 1; k <= kmax; ++k) risk_sets.push_back(std::vector<uint32_t>(ris_seq.begin(), ris_seq.begin() + k));
    auto ris_sel = evaluate_sets(d, risk_sets, Wsel, 555000 + (uint64_t)seed, materialize);
    auto quant = [&](const std::vector<uint32_t>& x, double q) { std::vector<double> v(x.begin(), x.end()); std::sort(v.begin(), v.end()); return v[(size_t)std::floor(q * (v.size() - 1))]; };
    double tlo = std::max(1.0, quant(ris_sel[0], 0.05)), thi = std::max(tlo * 1.01, quant(ris_sel[kmax - 1], 0.2));
    std::vector<double> taus; for (uint32_t t = 0; t < NT; ++t) taus.push_back(NT == 1 ? tlo : tlo * std::pow(thi / tlo, (double)t / (NT - 1)));
    size_t A = alphas.size(), NO = 1 + A + taus.size();
    std::vector<std::vector<uint32_t>> cur(NO);
    std::vector<std::vector<std::vector<uint32_t>>> chosen_sp(NO, std::vector<std::vector<uint32_t>>(kmax + 1));
    std::map<std::pair<size_t, long>, std::vector<uint32_t>> orc;
    Timer T2;
    for (uint32_t step = 1; step <= kmax; ++step) {
        std::vector<std::vector<uint32_t>> cand; std::vector<std::pair<size_t, uint32_t>> who;
        for (size_t o = 0; o < NO; ++o) for (uint32_t u : pool) {
            if (std::find(cur[o].begin(), cur[o].end(), u) != cur[o].end()) continue;
            std::vector<uint32_t> S = cur[o]; S.push_back(u); cand.push_back(S); who.push_back({o, u});
        }
        auto sp = evaluate_sets(d, cand, Wsel, 555000 + (uint64_t)seed, materialize);
        std::vector<double> best(NO, -1.0), bestm(NO, -1.0); std::vector<uint32_t> arg(NO, UINT32_MAX); std::vector<size_t> argc(NO, 0);
        std::vector<double> v(Wsel);
        for (size_t c = 0; c < cand.size(); ++c) {
            size_t o = who[c].first;
            for (uint32_t w = 0; w < Wsel; ++w) v[w] = sp[c][w];
            double mu = mean_of(v), sc = mu;
            if (o >= 1 && o <= A) { std::sort(v.begin(), v.end()); sc = cvar_sorted(v, alphas[o - 1]); }
            else if (o > A) { double tau = taus[o - 1 - A]; sc = 0; for (double x : v) sc += std::min(x, tau); sc /= Wsel; }
            if (sc > best[o] + 1e-12 || (std::fabs(sc - best[o]) <= 1e-12 && mu > bestm[o])) { best[o] = sc; bestm[o] = mu; arg[o] = who[c].second; argc[o] = c; }
        }
        for (size_t o = 0; o < NO; ++o) {
            cur[o].push_back(arg[o]); chosen_sp[o][step] = sp[argc[o]];
            if (std::find(ks.begin(), ks.end(), (long)step) != ks.end()) orc[{o, (long)step}] = cur[o];
        }
    }
    // Oracle-basis: for each (k, alpha) the best of all forward-simulated candidates and the RIS prefix by exact CVaR on Wsel
    std::map<std::pair<size_t, long>, std::vector<uint32_t>> obasis;
    for (long k : ks) for (size_t ai = 0; ai < A; ++ai) {
        double bv = -1; std::vector<uint32_t> bs;
        auto consider = [&](const std::vector<uint32_t>& x, const std::vector<uint32_t>& S) {
            std::vector<double> v(x.begin(), x.end()); std::sort(v.begin(), v.end()); double c = cvar_sorted(v, alphas[ai]);
            if (c > bv) { bv = c; bs = S; } };
        consider(ris_sel[k - 1], risk_sets[k - 1]);
        { std::vector<uint32_t> hs = hv.query((uint32_t)k, alphas[ai]).seeds; auto hsp = evaluate_sets(d, {hs}, Wsel, 555000 + (uint64_t)seed, materialize); consider(hsp[0], hs); }
        for (size_t o = 0; o < NO; ++o) consider(chosen_sp[o][k], std::vector<uint32_t>(cur[o].begin(), cur[o].begin() + k));
        obasis[{ai, k}] = bs;
    }
    double t_oracle = T2.sec();
    // evaluation on fresh worlds (disjoint stream)
    struct Rec { std::string method; long k; double alpha; std::vector<uint32_t> S; int h; double rho; };
    std::vector<Rec> recs;
    for (long k : ks) {
        std::vector<uint32_t> ris(ris_seq.begin(), ris_seq.begin() + k);
        recs.push_back({"RIS-IndRR", k, 0.0, ris, -1, NAN});
        recs.push_back({"Oracle-mean", k, 0.0, orc[{0, k}], -1, NAN});
        for (size_t ai = 0; ai < alphas.size(); ++ai) {
            auto ans = hv.query((uint32_t)k, alphas[ai]);
            recs.push_back({"RiskIndex-HV", k, alphas[ai], ans.seeds, (int)ans.h, ans.rho});
            if (with_fl) recs.push_back({fl_label, k, alphas[ai], fl.serves((uint32_t)k) ? fl.query((uint32_t)k, alphas[ai]).seeds : ans.seeds, -1, NAN});
            recs.push_back({"Oracle-CVaR", k, alphas[ai], orc[{ai + 1, k}], -1, NAN});
            recs.push_back({"Oracle-basis", k, alphas[ai], obasis[{ai, k}], -1, NAN});
        }
    }
    // k = 1 certificate: the 20 best-covered nodes outside the pool; no single node outside the pool can have a
    // CVaR above the largest of their means if the coverage ranking is right (CVaR <= mean).
    for (uint32_t i = NS; i < std::min<uint32_t>(NS + 20, g.n); ++i) recs.push_back({"PoolEdge", 1, 0.0, {order[i]}, -1, NAN});
    std::vector<std::vector<uint32_t>> sets; for (auto& r : recs) sets.push_back(r.S);
    auto sp = evaluate_sets(d, sets, Weval, 777000 + (uint64_t)seed, materialize);
    std::vector<double> v(Weval);
    for (size_t i = 0; i < recs.size(); ++i) {
        for (uint32_t w = 0; w < Weval; ++w) v[w] = sp[i][w];
        double mu = mean_of(v); std::sort(v.begin(), v.end());
        fprintf(out, "%s,%s,%ld,%u,%u,%s,%ld,%g,%.6f", gname.c_str(), model_name(model), seed, M, b, recs[i].method.c_str(), recs[i].k, recs[i].alpha, mu);
        for (double al : alphas) fprintf(out, ",%.6f", cvar_sorted(v, al));
        fprintf(out, ",%d,%.4f,%s\n", recs[i].h, recs[i].rho, vec_str(recs[i].S).c_str());
    }
    fflush(out); fclose(out);
    fprintf(stderr, "[headroom %s %s seed=%ld b=%u] index %.1fs oracle %.1fs pool=%zu\n", gname.c_str(), model_name(model), seed, b, t_index, t_oracle, pool.size());
    return 0;
}

int main(int argc, char** argv) {
    if (argc < 2) { fprintf(stderr, "usage: riskindex <oracle|fidelity|quality|scale|info> key=value ...\n"); return 1; }
    Args a; for (int i = 2; i < argc; ++i) { std::string s = argv[i]; auto p = s.find('='); if (p != std::string::npos) a.kv[s.substr(0, p)] = s.substr(p + 1); }
    std::string cmd = argv[1];
    try {
        if (cmd == "oracle") return cmd_oracle(a);
    if (cmd == "selftest") return cmd_selftest(a);
        if (cmd == "fidelity") return cmd_fidelity(a);
        if (cmd == "quality") return cmd_quality(a);
        if (cmd == "scale") return cmd_scale(a);
        if (cmd == "info") return cmd_info(a);
        if (cmd == "dist") return cmd_dist(a);
        if (cmd == "pair") return cmd_pair(a);
        if (cmd == "headroom") return cmd_headroom(a);
    } catch (std::exception& e) { fprintf(stderr, "error: %s\n", e.what()); return 2; }
    fprintf(stderr, "unknown command\n"); return 1;
}
