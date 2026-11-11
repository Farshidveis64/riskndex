# RiskIndex

Code for the ICDE 2027 submission "RiskIndex: Fast Tail-Risk Influence Maximization via Reusable Graph Indexing".

RiskIndex builds an index once per graph and then answers tail-risk influence queries: for a budget k and a
risk level alpha it returns a seed set chosen for the lower-tail CVaR of its spread, together with an
out-of-sample estimate of that CVaR, with a few array reads per query.

## Requirements

- C++17 compiler (tested with GCC 13.3); no external C++ libraries
- Python 3 (tested with 3.13) and `pip install -r requirements.txt`, for the analysis and figures

## Quick start

```bash
scripts/build.sh                  # builds bin/riskindex
bin/riskindex selftest            # correctness checks; prints "0 failures" and exits 0
```

## Reproduce the paper

```bash
scripts/download_data.sh          # SNAP graphs -> datasets/raw
scripts/run_all.sh                # all experiments (about 60 h on one core)
scripts/make_figures.sh           # tables and figures from results/
```

With the released raw results in `results/`, `scripts/ad_analyze.sh` recomputes every table and number of the
paper in about 30 seconds.

## Repository layout

| Path | Content |
|---|---|
| `src/` | C++17 implementation (header-only library and the driver `main.cpp`) |
| `scripts/` | build, data download and one script per experiment |
| `experiments/` | analysis and figure scripts (Python, TikZ) |
| `configs/parameters.md` | every parameter of every experiment |
| `AUDIT.md` | each claim of the paper -> result file -> script |
| `docs/DETAILS.md` | paper-to-code map, driver options, result files, method labels |

## Citation

```bibtex
@inproceedings{veisi2027riskindex,
  author    = {Farshid Veisi and Kamala Berahmand and Sona Taheri and Mahnaz Moradi and Mahdi Jalili},
  title     = {{RiskIndex}: Fast Tail-Risk Influence Maximization via Reusable Graph Indexing},
  booktitle = {Submitted to the 43rd IEEE International Conference on Data Engineering (ICDE)},
  year      = {2027}
}
```

## License

MIT (see `LICENSE`). The SNAP graphs are not redistributed; `scripts/download_data.sh` downloads them from
https://snap.stanford.edu/data.
