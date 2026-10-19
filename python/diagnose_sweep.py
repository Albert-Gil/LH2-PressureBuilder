#!/usr/bin/env python3
"""Independent equation-level audit of unsolved saved sweep points.

This diagnostic does not execute MATLAB or alter its accepted domain. Dense
sign-change scans remain finite numerical evidence, not a root-existence proof.
"""
from collections import Counter
from functools import lru_cache
from pathlib import Path
import argparse
import hashlib
import json

import CoolProp
from CoolProp.CoolProp import PropsSI
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Patch
import numpy as np
from scipy.io import loadmat
from scipy.optimize import brentq


ROOT = Path(__file__).resolve().parents[1]


class EquationAudit:
    def __init__(self, p):
        self.p = p
        self.cp = lambda output, key, value: PropsSI(output, 'P', p['P'], key, value, p['fluid'])
        self.hL = self.cp('H', 'T', p['T_L'])
        self.rhoL = self.cp('D', 'T', p['T_L'])
        self.muL = self.cp('V', 'T', p['T_L'])
        self.hf, self.hg = self.cp('H', 'Q', 0), self.cp('H', 'Q', 1)
        self.rf, self.rg = self.cp('D', 'Q', 0), self.cp('D', 'Q', 1)
        self.muf, self.mug = self.cp('V', 'Q', 0), self.cp('V', 'Q', 1)
        self.ts = self.cp('T', 'Q', 0)
        self.eos_limit = self.cp('Tmax', 'T', p['T_L'])

    def friction(self, re):
        lam = 64 / re
        turbulent = (-1.8 * np.log10((self.p['roughness'] / self.p['D'] / 3.7)**1.11 + 6.9/re))**-2
        s = np.clip((re - 2300)/1700, 0, 1)
        weight = s*s*(3-2*s)
        return (1-weight)*lam + weight*turbulent

    def state(self, flow, power):
        flow = np.atleast_1d(flow).astype(float)
        h = self.hL + power/flow
        rho, mu, temp = np.empty_like(h), np.empty_like(h), np.empty_like(h)
        mix = (h >= self.hf) & (h <= self.hg)
        x = (h[mix] - self.hf)/(self.hg-self.hf)
        rho[mix] = 1/((1-x)/self.rf + x/self.rg)
        mu[mix] = 1/((1-x)/self.muf + x/self.mug)
        temp[mix] = self.ts
        if np.any(~mix):
            props = np.asarray(PropsSI(['D', 'V', 'T'], 'P', self.p['P'],
                                       'H', h[~mix], self.p['fluid'])).reshape(-1, 3)
            rho[~mix], mu[~mix], temp[~mix] = props.T
        p = self.p
        area = np.pi*p['D']**2/4
        flux2 = (flow/area)**2
        fl = self.friction(flow*p['D']/area/self.muL)
        fg = self.friction(flow*p['D']/area/mu)
        friction = flux2/2*(fl*p['L_liq']/p['D']/self.rhoL + fg*p['L_vap']/p['D']/rho)
        minor = flux2/2*(p['K_L']/self.rhoL + p['K_g']/rho)
        acceleration = flux2*(1/rho - 1/self.rhoL)
        residual_at_zero_height = -p['g']*rho*p['H_elev'] - friction - minor - acceleration
        return residual_at_zero_height, temp

    @lru_cache(maxsize=256)
    def scan_grid(self, power, cap, samples):
        low = max(self.p['mdot_min'], power/(self.cp('H', 'T', cap)-self.hL)*(1+1e-9))
        high = self.p['mdot_hard_max']
        if low >= high:
            return np.array([]), np.array([])
        phase = power / (np.array([self.hf, self.hg])-self.hL)
        grid = np.unique(np.r_[np.geomspace(low, high, samples), phase[(phase > low) & (phase < high)]])
        residual, _ = self.state(grid, power)
        return grid, residual

    def roots(self, height, power, cap, samples):
        grid, base = self.scan_grid(float(power), float(cap), samples)
        values = base + self.p['g']*self.rhoL*height
        roots = list(grid[values == 0])
        for k in np.flatnonzero(values[:-1]*values[1:] < 0):
            root = brentq(lambda m: self.state(m, power)[0][0] + self.p['g']*self.rhoL*height,
                          grid[k], grid[k+1], xtol=1e-13)
            if not roots or min(abs(root-r) for r in roots) > 1e-9:
                roots.append(root)
        return sorted(roots)


def diagnose(directory, samples=512):
    directory = Path(directory)
    source = directory/'sweep_results.mat'
    s = loadmat(source, simplify_cells=True)
    audit = EquationAudit(s['p'])
    cap = min(float(s['p']['T_out_max']), audit.eos_limit)
    ok = np.asarray(s['ok'], bool)
    categories = np.full(ok.shape, 'accepted', dtype=object)
    records = []
    # Recheck an accepted point to catch unit/closure errors in this audit.
    control = s['sp']
    residue, temp = audit.state(control['mdot'], s['p']['Qdot_ref'])
    residual = float(residue[0] + s['p']['g']*audit.rhoL*s['p']['H_ref'])
    if abs(residual) > 1e-3 or abs(temp[0]-control['T_out']) > 1e-6:
        raise ValueError('Independent equations do not reproduce the saved reference point.')
    for row, col in np.argwhere(~ok):
        height, power = float(s['H_list'][row]), float(s['Q_list'][col])
        current = audit.roots(height, power, cap, samples)
        relaxed = audit.roots(height, power, audit.eos_limit, samples)
        relaxed_temperatures = [float(audit.state(m, power)[1][0]) for m in relaxed]
        # A root recovered by either grid inside the original domain is a
        # sampling discrepancy, not evidence for temperature-cap exclusion.
        if not current:
            current = [m for m, t in zip(relaxed, relaxed_temperatures) if t <= cap+1e-6]
        if str(s['status'][row, col]) != 'no_root_in_domain':
            category = 'other_recorded_failure'
        elif current:
            category = 'root_found_on_dense_scan'
        elif relaxed:
            category = 'temperature_cap_excluded'
        else:
            category = 'no_root_found_to_eos_limit'
        categories[row, col] = category
        records.append({'row': int(row), 'column': int(col), 'height_m': height,
                        'power_W': power, 'recorded_status': str(s['status'][row, col]),
                        'diagnosis': category, 'roots_at_configured_cap_kg_s': current,
                        'roots_at_eos_limit_kg_s': relaxed,
                        'relaxed_root_temperatures_K': relaxed_temperatures})
    counts = dict(Counter(categories.ravel()))
    result = {'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
              'coolprop_version': CoolProp.__version__, 'samples_per_domain': samples,
              'configured_temperature_cap_K': cap, 'eos_temperature_limit_K': audit.eos_limit,
              'flow_hard_max_kg_s': float(s['p']['mdot_hard_max']),
              'control_pressure_residual_Pa': residual, 'counts': counts,
              'caveat': 'Independent Python equation audit, not MATLAB execution. Finite sign-change scans are not a proof of no root.',
              'records': records}
    (directory/'sweep_diagnostics.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    plot_diagnosis(s, categories, directory)
    write_report(result, directory)
    return result


def plot_diagnosis(s, categories, directory):
    names = ['accepted', 'temperature_cap_excluded', 'no_root_found_to_eos_limit',
             'root_found_on_dense_scan', 'other_recorded_failure']
    palette = ['#087f8c', '#e9b44c', '#b94e48', '#7b61a8', '#718096']
    present = [k for k, name in enumerate(names) if np.any(categories == name)]
    codes = np.zeros(categories.shape, int)
    for code, k in enumerate(present):
        codes[categories == names[k]] = code
    fig, ax = plt.subplots(figsize=(11, 7), constrained_layout=True)
    ax.pcolormesh(s['Q_list']/1000, s['H_list'], codes, shading='nearest',
                  cmap=ListedColormap([palette[k] for k in present]),
                  norm=BoundaryNorm(np.arange(len(present)+1)-.5, len(present)))
    labels = {'accepted': 'Accepted in saved MATLAB sweep',
              'temperature_cap_excluded': 'Root found only above configured temperature cap',
              'no_root_found_to_eos_limit': 'No root found up to EOS temperature limit',
              'root_found_on_dense_scan': 'Root found by denser in-domain scan',
              'other_recorded_failure': 'Other recorded failure'}
    handles = [Patch(facecolor=palette[k], label=f'{labels[names[k]]} ({np.count_nonzero(categories == names[k])})') for k in present]
    fig.legend(handles=handles, loc='outside lower center', fontsize=11)
    ax.set(xlabel='Absorbed heater power [kW]', ylabel='Liquid height [m]',
           title='Unsolved-region diagnosis\nIndependent equation audit; finite root scans')
    folder = directory/'figures'
    folder.mkdir(parents=True, exist_ok=True)
    fig.savefig(folder/'sweep_failure_diagnosis_H_Q.png', dpi=180)
    plt.close(fig)


def write_report(result, directory):
    lines = ['# Sweep domain diagnosis', '',
             result['caveat'], '',
             f"Source MAT SHA-256: `{result['source_sha256']}`.", '',
             f"CoolProp {result['coolprop_version']}; {result['samples_per_domain']} logarithmic samples per domain, plus phase boundaries. Flow searched to {result['flow_hard_max_kg_s']:g} kg/s. Configured cap: {result['configured_temperature_cap_K']:g} K; relaxed cap: EOS limit {result['eos_temperature_limit_K']:g} K.", '',
             '| Category | Points |', '|---|---:|']
    lines += [f'| `{name}` | {count} |' for name, count in result['counts'].items()]
    lines += ['', '## Representative unsolved points', '',
              '| H [m] | Q [kW] | Diagnosis | Relaxed roots [kg/min] | Outlet T [K] |',
              '|---:|---:|---|---|---|']
    for name in sorted(set(r['diagnosis'] for r in result['records'])):
        cases = [r for r in result['records'] if r['diagnosis'] == name]
        for idx in sorted({0, len(cases)//2, len(cases)-1}):
            r = cases[idx]
            flows = ', '.join(f'{60*m:.6f}' for m in r['roots_at_eos_limit_kg_s']) or 'None found'
            temps = ', '.join(f'{t:.3f}' for t in r['relaxed_root_temperatures_K']) or 'Undefined'
            lines.append(f"| {r['height_m']:.6f} | {r['power_W']/1000:.6f} | `{name}` | {flows} | {temps} |")
    lines += ['', '## Interpretation', '',
              '- `temperature_cap_excluded`: the denser scan finds no configured-domain root but finds a pressure-balanced root when the cap is relaxed. These out-of-domain states are diagnostic only; they are not accepted hardware operating points.',
              '- `no_root_found_to_eos_limit`: neither scan found a sign-changing root. This does not establish physical stagnation; roots outside the EOS/flow bounds or missed by sampling are not ruled out.',
              '- `root_found_on_dense_scan`: investigate sampling/bounds in the original solver before interpreting the saved failure.',
              '- Source data, configured cap, and production solver remain unchanged.', '']
    (directory/'sweep_diagnostics.md').write_text('\n'.join(lines), encoding='utf-8')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory', type=Path, default=ROOT/'results')
    parser.add_argument('--samples', type=int, default=512)
    args = parser.parse_args()
    if args.samples < 48:
        parser.error('--samples must be at least 48')
    summary = diagnose(args.directory, args.samples)
    print(json.dumps(summary['counts'], indent=2))
