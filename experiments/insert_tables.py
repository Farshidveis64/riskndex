#!/usr/bin/env python3
"""Insert the generated table bodies (results/tables_ad/*.tex) into the paper (main.tex): each body
replaces the rows between the header's \\midrule and \\bottomrule of the table identified by its header line.
Usage: insert_tables.py <main.tex> [results_dir]"""
import sys, os
tex = sys.argv[1]; R = sys.argv[2] if len(sys.argv) > 2 else 'results'
s = open(tex).read()
specs = [('quality.tex', 'Graph & Model & $\\E[\\sigma]$ & $\\CVaR_{0.05}$'),
         ('ablation.tex', 'Variant & G & B & S & M & V & R & O & E'),
         ('systems.tex', 'Type & Graph & $b$ & sketch & index &'),
         ('theta.tex', '$\\theta$ & Tier A (\\%) & Tier B gain (\\%) & layer (\\%)')]
for fn, head in specs:
    f = f'{R}/{os.environ.get("TABDIR", "tables_ad")}/{fn}'
    if not os.path.exists(f) or head not in s: print('skip', fn); continue
    body = open(f).read().rstrip('\n') + '\n'
    i = s.index(head); j = s.index('\\midrule\n', i) + len('\\midrule\n'); k = s.index('\\bottomrule', j)
    s = s[:j] + body + s[k:]
    print('inserted', fn, body.count('\\\\'), 'rows')
open(tex, 'w').write(s)
