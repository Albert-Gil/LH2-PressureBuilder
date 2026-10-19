#!/usr/bin/env python3
"""Batch plots for the MAT files in results/ (optional; MATLAB plots the same figures)."""
from pathlib import Path
import argparse
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Patch
import numpy as np
from scipy.io import loadmat

RESULTS = Path(__file__).resolve().parent.parent / "results"
plt.rcParams.update({"font.size": 18, "axes.titlesize": 19, "legend.fontsize": 14})


def load_results(path):
    data = loadmat(path, simplify_cells=True)
    if "ok" not in data or "metadata" not in data:
        raise ValueError(f"{path.name} predates the repaired model; rerun its MATLAB entry point.")
    failed = np.size(data["ok"]) - np.count_nonzero(data["ok"])
    if failed:
        warnings.warn(f"{path.name}: {failed} failed points are shown as gaps.", stacklevel=2)
    return data


def save(fig, directory, name):
    folder = Path(directory) / "figures"
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / f"{name}.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def decorate(ax, xlabel, ylabel, title=""):
    ax.set(xlabel=xlabel, ylabel=ylabel, title=title)
    ax.grid(True, alpha=0.25)


def plot_scenario(kind, directory=RESULTS):
    if kind not in ("SOC", "constH"):
        raise ValueError("kind must be SOC or constH")
    directory = Path(directory)
    s = load_results(directory / f"{kind}_scenario_results.mat")
    t = np.atleast_1d(s["tmin"])
    h = np.atleast_1d(s["H"])
    labels = [str(x) for x in np.atleast_1d(s["labels"]).ravel()]
    ncase = len(labels)
    shape = (len(t), ncase)
    q = np.asarray(s["Qdot"]).reshape(shape) / 1000
    ok = np.asarray(s["ok"], dtype=bool).reshape(shape)
    mdot = np.where(ok, np.asarray(s["mdot"]).reshape(shape), np.nan) * 60
    tout = np.where(ok, np.asarray(s["T_out"]).reshape(shape), np.nan)
    colors = plt.get_cmap("tab10")
    context = (f"prescribed volume fill {100*float(s['SOC0']):g}% to {100*float(s['SOC1']):g}%"
               if kind == "SOC" else f"constant height {float(s['H_const']):.2f} m")
    fig, axes = plt.subplots(2, 1, figsize=(13, 8), constrained_layout=True)
    fig.suptitle(f"LH2 vaporizer - {context}", fontweight="bold")
    if kind == "SOC":
        axes[0].plot(t, 100*np.atleast_1d(s["SOC"]), color="black")
        right = axes[0].twinx()
        right.plot(t, h, color="C0")
        right.set_ylabel("Liquid height [m]")
        decorate(axes[0], "Time [min]", "Volume fill [%]", "Prescribed tank level")
    else:
        axes[0].plot(t, h, color="black")
        decorate(axes[0], "Time [min]", "Liquid height [m]", "Prescribed tank level")
    for k, label in enumerate(labels):
        axes[1].plot(t, q[:, k], color=colors(k % 10), label=label)
    decorate(axes[1], "Time [min]", "Absorbed heater power [kW]", "Shared power schedules")
    axes[1].margins(y=0.08)  # Autoscaling includes every case, including constant power.
    fig.legend(*axes[1].get_legend_handles_labels(), loc="outside right center")
    outputs = [save(fig, directory, f"{kind}_scenario_drivers")]
    fig, axes = plt.subplots(2, 1, figsize=(13, 8), constrained_layout=True)
    fig.suptitle(f"LH2 vaporizer - {context}", fontweight="bold")
    for k, label in enumerate(labels):
        axes[0].plot(t, mdot[:, k], color=colors(k % 10), label=label)
        axes[1].plot(t, tout[:, k], color=colors(k % 10), label=label)
    decorate(axes[0], "Time [min]", "Mass flow [kg/min]", "Circulation through the heater")
    decorate(axes[1], "Time [min]", "Outlet temperature [K]", "Heater outlet / adiabatic return")
    fig.legend(*axes[0].get_legend_handles_labels(), loc="outside right center")
    outputs.append(save(fig, directory, f"{kind}_scenario_mdot_Tout"))
    return outputs


def plot_sweep(directory=RESULTS):
    directory = Path(directory)
    s = load_results(directory / "sweep_results.mat")
    h, q = np.atleast_1d(s["H_list"]), np.atleast_1d(s["Q_list"])/1000
    shape = (len(h), len(q))
    ok = np.asarray(s["ok"], dtype=bool).reshape(shape)
    m = np.where(ok, np.asarray(s["mdot_map"]).reshape(shape), np.nan)*60
    temp = np.where(ok, np.asarray(s["Tout_map"]).reshape(shape), np.nan)
    x = np.where(ok, np.asarray(s["x_map"]).reshape(shape), np.nan)
    ih, iq = int(s["iHref"])-1, int(s["iQref"])-1
    outputs = []
    for data, name, label in [(m, "mdot_map_H_Q", "Mass flow [kg/min]"),
                              (temp, "T_out_map_H_Q", "Outlet temperature [K]")]:
        fig, ax = plt.subplots(figsize=(8.5, 6), constrained_layout=True)
        ax.set_facecolor('#e5e7eb')
        # pcolormesh also supports one-row or one-column custom sweeps.
        mesh = ax.pcolormesh(q, h, np.ma.masked_invalid(data), shading="nearest", cmap="viridis")
        fig.colorbar(mesh, ax=ax, label=label)
        ax.set(xlabel="Absorbed heater power [kW]", ylabel="Liquid height [m]",
               title=f"Horizontal heater / vertical return - {label}")
        if np.any(~ok):
            ax.legend(handles=[Patch(facecolor='#e5e7eb', label='No accepted solution')],
                      loc='lower right', fontsize=10)
        outputs.append(save(fig, directory, name))
    specs = [(q, m[ih], "mdot_vs_Q", "Absorbed heater power [kW]", "Mass flow [kg/min]"),
             (h, m[:, iq], "mdot_vs_H", "Liquid height [m]", "Mass flow [kg/min]"),
             (q, x[ih], "x_out_vs_Q", "Absorbed heater power [kW]", "Outlet vapor fraction"),
             (q, temp[ih], "T_out_vs_Q", "Absorbed heater power [kW]", "Outlet temperature [K]")]
    for xs, ys, name, xlabel, ylabel in specs:
        fig, ax = plt.subplots(figsize=(8.5, 5.5), constrained_layout=True)
        ax.plot(xs, ys, "-o", ms=3, label="Circulation solution")
        title = f"Liquid height = {h[ih]:.3f} m"
        if name == "mdot_vs_Q":
            ax.plot(q, np.atleast_1d(s["ceiling_vs_Q"])*60, "--",
                    label="Complete-vaporization boundary Q/(hg-hL)")
            ax.legend()
        elif name == "mdot_vs_H":
            title = f"Absorbed heater power = {q[iq]:.3g} kW"
        elif name == "T_out_vs_Q":
            ax.axhline(s["sp"]["T_sat"], ls="--", color="gray", label="Saturation temperature")
            ax.legend()
        elif name == "x_out_vs_Q":
            ax.set_ylim(0, 1.05)
        decorate(ax, xlabel, ylabel, title)
        outputs.append(save(fig, directory, name))
    outputs.append(plot_sweep_validity(s, directory))
    return outputs


def plot_sweep_validity(s, directory=RESULTS):
    """Render saved solver categories; never infer failure causes from NaNs."""
    h, q = np.atleast_1d(s['H_list']), np.atleast_1d(s['Q_list']) / 1000
    status = np.asarray(s['status']).reshape(len(h), len(q))
    categories = sorted(set(status.ravel()), key=lambda x: (x != 'solved', x))
    codes = np.zeros(status.shape, dtype=int)
    palette = ['#087f8c', '#e9b44c', '#b94e48', '#7b61a8', '#718096', '#d95f02']
    for k, category in enumerate(categories):
        codes[status == category] = k
    colors = [palette[k % len(palette)] for k in range(len(categories))]
    fig, ax = plt.subplots(figsize=(10, 6.5), constrained_layout=True)
    ax.pcolormesh(q, h, codes, shading='nearest', cmap=ListedColormap(colors),
                  norm=BoundaryNorm(np.arange(len(categories) + 1) - .5, len(categories)))
    handles = [Patch(facecolor=colors[k], label=f'{name}: {np.count_nonzero(status == name):,}')
               for k, name in enumerate(categories)]
    fig.legend(handles=handles, loc='outside lower center', ncol=min(2, len(handles)), fontsize=12)
    ax.set(xlabel='Absorbed heater power [kW]', ylabel='Liquid height [m]',
           title='Sweep validity: recorded solver status\nUnsolved does not mean zero flow')
    return save(fig, directory, 'sweep_validity_H_Q')


def plot_power_decay(directory=RESULTS):
    directory = Path(directory)
    s = load_results(directory / "power_decay_results.mat")
    t = np.atleast_1d(s["t"])/60
    ok = np.atleast_1d(s["ok"]).astype(bool)
    fig, axes = plt.subplots(4, 1, figsize=(9, 10), constrained_layout=True)
    fig.suptitle("Quasi-steady power decay - imposed pressure and inlet temperature", fontweight="bold")
    axes[0].plot(t, np.atleast_1d(s["Qdot"])/1000)
    right = axes[0].twinx()
    right.plot(t, np.atleast_1d(s["H"]), color="C1")
    right.set_ylabel("Liquid height [m]")
    decorate(axes[0], "Time [min]", "Absorbed heater power [kW]")
    for ax, field, scale, label in [(axes[1], "mdot", 60, "Mass flow [kg/min]"),
                                    (axes[2], "T_out", 1, "Outlet temperature [K]"),
                                    (axes[3], "x_out", 1, "Outlet vapor fraction")]:
        ax.plot(t, np.where(ok, np.atleast_1d(s[field]), np.nan)*scale)
        decorate(ax, "Time [min]", label)
    axes[2].axhline(s["T_sat"], ls="--", color="gray")
    axes[3].set_ylim(0, 1.05)
    return [save(fig, directory, "power_decay_transient")]


def plot_all(directory=RESULTS):
    return (plot_sweep(directory) + plot_power_decay(directory)
            + plot_scenario("SOC", directory) + plot_scenario("constH", directory))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, default=RESULTS)
    args = parser.parse_args()
    for output in plot_all(args.directory):
        print(output)
