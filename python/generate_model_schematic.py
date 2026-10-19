#!/usr/bin/env python3
"""Generate an exact vector drawing of the model topology (not to scale).

Annotations use the parameter struct in the tracked reference sweep dataset.
"""
from pathlib import Path
import argparse

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle, FancyArrowPatch, FancyBboxPatch
from scipy.io import loadmat


ROOT = Path(__file__).resolve().parents[1]


def generate(directory=ROOT/'docs'):
    p = loadmat(ROOT/'results'/'sweep_results.mat', simplify_cells=True)['p']
    out = Path(directory)
    out.mkdir(parents=True, exist_ok=True)
    navy, blue, red, muted = '#17324d', '#087f8c', '#bd4f36', '#566879'
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 12,
                         'svg.fonttype': 'none'})
    fig, ax = plt.subplots(figsize=(17, 8), facecolor='white')
    ax.set(xlim=(-.5, 20), ylim=(-.3, 9.2), aspect='equal')
    ax.axis('off')
    ax.text(.2, 8.85, 'LH₂ PRESSURE-BUILDER LOOP', fontsize=23, weight='bold', color=navy)
    ax.text(.2, 8.28, 'Model schematic • circular tank cross-section • not to scale', fontsize=12, color=muted)

    cx, cy, radius = 4.5, 5, 2.5
    bottom, top = cy-radius, cy+radius
    level = bottom + 2*radius*p['H_ref']/p['D_tank']
    shell = Circle((cx, cy), radius, facecolor='#f4f8fc', edgecolor='none')
    ax.add_patch(shell)
    liquid = Rectangle((cx-radius, bottom), 2*radius, level-bottom, facecolor='#c6e8ec', edgecolor='none')
    liquid.set_clip_path(shell)
    ax.add_patch(liquid)
    surface, = ax.plot([cx-radius, cx+radius], [level, level], color=blue, lw=1.7)
    surface.set_clip_path(shell)
    ax.add_patch(Circle((cx, cy), radius, facecolor='none', edgecolor=navy, lw=2.5))
    ax.text(cx, 6.65, 'ULLAGE', ha='center', weight='bold', color=navy)
    ax.text(cx, 6.15, f'P = {p["P"]/1e5:.3f} bar abs.', ha='center', color=navy)
    ax.text(cx, 4.5, 'LIQUID HYDROGEN', ha='center', weight='bold', color=navy)
    ax.text(cx, 3.97, rf'$T_L = {p["T_L"]:g}$ K  |  $h_L,\,\rho_L$', ha='center', color=navy)
    ax.text(cx, 3.4, 'Horizontal cylindrical tank', ha='center', fontsize=10, color=muted)

    # The heater and takeoff share a datum, matching the implemented head.
    ax.plot([cx, 7.1], [bottom, bottom], color=blue, lw=4, solid_capstyle='round')
    ax.plot([9.1, 10.4, 10.4, cx, cx], [bottom, bottom, top, top, 7.02],
            color=red, lw=4, solid_capstyle='round', solid_joinstyle='round')
    heater = FancyBboxPatch((7.1, bottom-.4), 2, .8, boxstyle='round,pad=0.05,rounding_size=.08',
                           facecolor='#fff0e9', edgecolor=red, lw=2)
    ax.add_patch(heater)
    for x in [7.35, 7.65, 7.95, 8.25, 8.55, 8.85]:
        ax.plot([x, x], [bottom-.31, bottom+.31], color=red, lw=1)
    ax.text(8.1, 3.32, 'HORIZONTAL HEATER', ha='center', fontsize=11, weight='bold', color=red)
    ax.text(8.1, 1.54, 'Net absorbed heat', ha='center', fontsize=11, color=red)
    ax.annotate('', xy=(8.1, 2.04), xytext=(8.1, 1.73),
                arrowprops={'arrowstyle': '-|>', 'color': red, 'lw': 2})
    ax.text(8.1, .88, rf'$\dot Q_{{ref}} = {p["Qdot_ref"]/1000:g}$ kW', ha='center', color=red, fontsize=13, weight='bold')

    def arrow(start, end, color):
        ax.add_patch(FancyArrowPatch(start, end, arrowstyle='-|>', mutation_scale=19,
                                     color=color, lw=2.2, zorder=5))
    arrow((5.5, bottom), (6.35, bottom), blue)
    arrow((10.4, 3.65), (10.4, 4.6), red)
    arrow((8.35, top), (7.45, top), red)
    arrow((cx, 7.35), (cx, 7.02), red)
    ax.text(5.8, 1.82, 'Liquid supply', ha='center', color=blue, fontsize=11)
    ax.text(9.85, 5.3, 'Adiabatic return', rotation=90, ha='center', color=red)
    ax.text(8.4, 7.82, r'$\dot m,\ h_{out},\ T_{out},\ x_{out}$', ha='center', fontsize=12, color=red)

    def dimension(x, y0, y1, label, color=muted):
        ax.annotate('', xy=(x, y1), xytext=(x, y0),
                    arrowprops={'arrowstyle': '<->', 'color': color, 'lw': 1.2})
        for y in (y0, y1):
            ax.plot([x-.13, x+.13], [y, y], color=color, lw=1)
        ax.text(x-.23, (y0+y1)/2, label, rotation=90, va='center', ha='center', color=color, fontsize=11)
    dimension(.35, bottom, top, rf'$D_{{tank}} = {p["D_tank"]:g}$ m')
    dimension(1.35, bottom, level, rf'$H_{{ref}} = {p["H_ref"]:g}$ m', blue)
    dimension(11.85, bottom, top, rf'$H_{{elev}} = {p["H_elev"]:g}$ m')
    ax.plot([.15, 12.1], [bottom, bottom], color=muted, ls='--', lw=.7, zorder=0)
    ax.text(.3, 1.15, 'Takeoff / heater elevation datum', fontsize=10, color=muted)

    card = FancyBboxPatch((13.1, 2.28), 6.25, 5.35, boxstyle='round,pad=.2,rounding_size=.12',
                          facecolor='#f4f8fc', edgecolor='#d4dfe9', lw=1)
    ax.add_patch(card)
    ax.text(13.4, 7.1, 'REFERENCE MODEL INPUTS', weight='bold', fontsize=13, color=navy)
    rows = [f'Pipe ID, D       {p["D"]*1000:.2f} mm',
            f'Supply, Lliq     {p["L_liq"]:g} m',
            f'Return, Lvap     {p["L_vap"]:g} m',
            f'Losses, KL / Kg  {p["K_L"]:g} / {p["K_g"]:g}',
            f'Capacity         {p["vol_L"]/1000:g} m³',
            f'Cylinder length  {p["tank_length"]:.3f} m']
    for i, line in enumerate(rows):
        ax.text(13.4, 6.49-i*.49, line, fontfamily='DejaVu Sans Mono', fontsize=11, color=navy)
    ax.text(13.4, 2.72, 'Effective losses include unresolved\nheater passages and fittings.', color=muted, fontsize=10, linespacing=1.5)
    ax.text(.2, .02, 'Prescribed P and TL • uniform outlet state in return • no tank-pressure dynamics', color=muted, fontsize=11)
    fig.subplots_adjust(left=.02, right=.99, bottom=.04, top=.98)
    paths = [out/'model_schematic.svg', out/'model_schematic.png']
    for path in paths:
        fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig)
    return paths


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory', type=Path, default=ROOT/'docs')
    for path in generate(parser.parse_args().directory):
        print(path)
