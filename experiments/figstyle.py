"""Shared matplotlib style for all paper figures (IEEE two-column): muted palette, dashed grids,
sans-serif plot text, serif subcaptions below the panels.

Palette: see C below (validated); every figure must remain
readable in grayscale, so series are also distinguished by marker/linestyle.
"""
import matplotlib, matplotlib.font_manager as fm
matplotlib.use('Agg')
import matplotlib.pyplot as plt

for p in ['/usr/share/texmf/fonts/opentype/public/tex-gyre/texgyretermes-regular.otf',
          '/usr/share/texmf/fonts/opentype/public/tex-gyre/texgyretermes-bold.otf',
          '/usr/share/texmf/fonts/opentype/public/tex-gyre/texgyretermes-italic.otf',
          '/usr/share/texmf/fonts/opentype/public/tex-gyre/texgyretermes-bolditalic.otf']:
    try: fm.fontManager.addfont(p)
    except Exception: pass

COLW = 3.5    # IEEE column width (in)
TEXTW = 7.16  # IEEE text width (in)
ROWH = 0.86   # height of a 1x4 result row (in)

# Palette (validated with the dataviz validator, light surface):
#   three categorical hues blue/orange/aqua pass every gate all-pairs (CVD dE >= 9.2,
#   normal-vision dE >= 24).  More than three hues cannot pass all-pairs, so further
#   series are de-emphasized in gray and carried by marker + dash pattern (emphasis
#   encoding).  Every series also has a distinct marker/linestyle, so all figures
#   stay readable in grayscale print.
C = dict(purple='#7b5fb0', salmon='#ee9a6a', sage='#6e9e83', sky='#8db8dd', blue='#4a86c5', orange='#ee9a6a',
         green='#6e9e83', aqua='#6e9e83', red='#d9705a', yellow='#e3b04b', magenta='#d27ba6',
         black='#2b2b2b', ink2='#5b5b5b', gray='#9a9a9a', light='#e3e3e3', grid='#cfcfcf')
# Muted palette after the visual language of recent ML venue papers: "ours" in a distinct purple drawn
# on top, baselines in salmon / sage / sky / gray.  Every series also has its own marker and dash
# pattern, so the figures stay readable in grayscale and for colour-vision deficiencies.
# method -> (color, marker, linestyle)
METHOD = {
    'RiskIndex':        (C['purple'], 'P', '--'),
    'DirectCVaR-GRR':   (C['salmon'], 's', '-.'),
    'SG':               (C['sage'],   '^', ':'),
    'Mean-GRR':         (C['sky'],    'o', '-'),
    'RIS-IndRR':        (C['gray'],   'v', (0, (1, 1.2))),
    'DirectCVaR-IndRR': (C['gray'],   'D', (0, (1, 1.2))),
    'SG-cached':        (C['gray'],   'X', (0, (3, 1, 1, 1))),
    'Oracle-Mean':      (C['black'],  'x', (0, (4, 1.5))),
    'Oracle':           (C['black'],  '*', '-'),
    'GRR':              (C['purple'], 'P', '--'),
    'IndRR':            (C['salmon'], 's', '-.'),
}
# graphs (sensitivity / fidelity): color = graph, marker fill = model
GRAPHC = {'wiki': C['purple'], 'email': C['salmon'], 'epinions': C['sage']}
MODELC = {'WC': C['purple'], 'LT': C['salmon'], 'IC': C['sage']}
LABEL = {
    'RiskIndex': 'RiskIndex', 'DirectCVaR-GRR': 'Direct-CVaR (GRR)', 'DirectCVaR-IndRR': 'Direct-CVaR (ind. RR)',
    'Mean-GRR': 'RIS-Mean (GRR)', 'RIS-IndRR': 'RIS-Mean (ind. RR)', 'Oracle-Mean': 'Exhaustive mean-opt.', 'SG': 'SG [10]', 'SG-cached': 'SG-cached',
    'Oracle': 'Exhaustive CVaR-opt.', 'GRR': 'GRR (grouped)', 'IndRR': 'independent RR',
}

plt.rcParams.update({
    # sans-serif inside the plots (ticks, labels, legends), serif subcaptions under the panels (see sub())
    'font.family': 'sans-serif', 'font.sans-serif': ['DejaVu Sans', 'Liberation Sans', 'Arial'],
    'mathtext.fontset': 'dejavusans',
    'font.size': 6.5, 'axes.labelsize': 6.5, 'axes.titlesize': 6.5, 'legend.fontsize': 5.8, 'xtick.labelsize': 5.8, 'ytick.labelsize': 5.8,
    'axes.linewidth': 0.5, 'axes.edgecolor': '#333333', 'xtick.color': '#333333', 'ytick.color': '#333333',
    'xtick.major.width': 0.5, 'ytick.major.width': 0.5, 'xtick.major.size': 2.0, 'ytick.major.size': 2.0,
    'lines.linewidth': 1.0, 'lines.markersize': 3.2, 'legend.frameon': False, 'legend.handlelength': 2.4,
    'axes.spines.top': True, 'axes.spines.right': True, 'figure.dpi': 150, 'savefig.dpi': 300,
    'pdf.fonttype': 42, 'ps.fonttype': 42, 'axes.grid': False,
})
SERIF = dict(family='serif', fontsize=7.5, math_fontfamily='stix')  # subcaptions: match the paper's Times text
for p in ['/usr/share/texmf/fonts/opentype/public/tex-gyre/texgyretermes-regular.otf']:
    try: fm.fontManager.addfont(p)
    except Exception: pass
plt.rcParams['font.serif'] = ['TeX Gyre Termes', 'Liberation Serif', 'Times New Roman', 'DejaVu Serif']

def style(ax, grid='both'):
    """Light dashed grid on both axes, drawn under the data."""
    if grid: ax.grid(axis=grid, color=C['grid'], linewidth=0.4, linestyle=(0, (3, 2)), zorder=0)
    ax.set_axisbelow(True)

def sub(ax, text, y=-0.42):
    """Subcaption '(a) ...' centred below the panel, in the paper's serif font."""
    ax.text(0.5, y, text, transform=ax.transAxes, ha='center', va='top', **SERIF)

# Shared frame (inches, figure coordinates) of the 1x4 result rows Fig. 6 (fig_fidelity.py) and Fig. 7
# (fig_sens_v2.py): the union of their tight boxes plus 0.01 in of padding. Both rows use the same figure size and
# grid, so with this frame the two PDFs have identical size and identical panel positions, and LaTeX scales them
# by the same factor (same height and font size on the page).
ROW_BBOX = (0.5064, -0.3018, 6.5475, 1.0381)  # with ROWH = 0.86

def save(fig, name, figdir='figures', bbox=None):
    """Save PDF and PNG; bbox=None crops tightly, bbox=(x0, y0, x1, y1) in inches uses a fixed frame."""
    import os; os.makedirs(figdir, exist_ok=True)
    if os.environ.get('PRINT_BBOX'):
        bb = fig.get_tightbbox(fig.canvas.get_renderer()); print('TIGHTBBOX', name, [round(v, 4) for v in (bb.x0, bb.y0, bb.x1, bb.y1)])
    if bbox is None:
        kw = dict(bbox_inches='tight', pad_inches=0.01)
    else:
        from matplotlib.transforms import Bbox
        kw = dict(bbox_inches=Bbox.from_extents(*bbox))
    fig.savefig(f'{figdir}/{name}.pdf', **kw)
    fig.savefig(f'{figdir}/{name}.png', dpi=300, **kw)
    print('saved', name)
