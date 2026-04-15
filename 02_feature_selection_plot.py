# Written by Mohd Ibrahim
# Technical University of Munich
import argparse
import pickle
import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
from matplotlib import gridspec


def parse_args():
    parser = argparse.ArgumentParser(description="Plot ΔG, H-bonds, and ion interactions from MD data")
    
    # General inputs    
    parser.add_argument('-p', dest='pkl_file', type=str, default='output_features.pkl', help='input pkl file')
    parser.add_argument('-c', dest='gro_file', type=str, default='mol.gro', help='input gro file')
    parser.add_argument('-s', dest='tpr_file', type=str, default='npt.tpr', help='input tpr file')
    
    # Thresholds
    parser.add_argument('-acceptor_th', dest='acceptor_threshold', type=float, default=40, help='acceptor threshold')
    parser.add_argument('-donor_th', dest='donor_threshold', type=float, default=40, help='donor threshold')
    parser.add_argument('-dG_th', dest='dG_threshold', type=float, default=0.4, help='ΔG threshold')
    parser.add_argument('-ion_th', dest='ion_threshold', type=float, default=25, help='difference between number of cations and anions')

    # Figure options
    parser.add_argument('-o', dest='out_file', type=str, default="features.png", help="Output image file (e.g. plot.png). If None, shows figure.")
    parser.add_argument('-figsize', dest='figsize', type=float, nargs=2, default=(12, 12), help="Figure size (width height)")
    
    return parser.parse_args()


def main():
    
    args = parse_args()
    pkl_file           = args.pkl_file
    gro_file           = args.gro_file
    tpr_file           = args.tpr_file
    acceptor_threshold = args.acceptor_threshold
    donor_threshold    = args.donor_threshold
    dG_threshold       = args.dG_threshold
    ion_threshold      = args.ion_threshold
    figsize            = args.figsize
    out_file           = args.out_file
    
    # ---------------- Load data ---------------- #
    with open(pkl_file, "rb") as f:
        data = pickle.load(f)

    # ---------------- Universe ---------------- #
    u = mda.Universe(tpr_file, gro_file)
    
    frag1, frag2 = u.atoms.fragments[0], u.atoms.fragments[1]
    Resid_frag1, Resid_frag2 = frag1.resids, frag2.resids
    Index_frag1, Index_frag2 = frag1.indices, frag2.indices

    # ---------------- Figure Setup ---------------- #
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(nrows=4, ncols=1, hspace=0.8)
    
    ax1 = fig.add_subplot(gs[0, 0])  # ΔG panel
    ax2 = fig.add_subplot(gs[1, 0])  # H-bond acceptors
    ax3 = fig.add_subplot(gs[2, 0])  # H-bond donors
    ax4 = fig.add_subplot(gs[3, 0])  # Ion binding

    # =====================================================
    # Panel 1: ΔG
    # =====================================================
    dG1, std1, resname1 = [], [], []
    dG2, std2, resname2 = [], [], []

    for key, values in data["sasa_dG"].items():
        resid = int(key)
        if resid in Resid_frag1:
            dG1.append(values['dG'][0])
            std1.append(values['dG'][1])
            resname1.append(values['resname'])
        else:
            dG2.append(values['dG'][0])
            std2.append(values['dG'][1])
            resname2.append(values['resname'])

    # Safety: if one chain has no data, still handle plotting
    x1 = np.arange(len(dG1))
    x2 = np.arange(len(dG2)) + (x1[-1] + 1 if len(x1) > 0 else 0)

    if len(dG1) > 0:
        ax1.bar(x1, dG1, yerr=std1, label="Chain A", edgecolor="k", lw=1, capsize=5)
    if len(dG2) > 0:
        ax1.bar(x2, dG2, yerr=std2, label="Chain B", edgecolor="k", lw=1, capsize=5, fc='tab:red')

    all_names = resname1 + resname2
    all_x = np.concatenate([x1, x2]) if len(all_names) > 0 else []
    if len(all_x) > 0:
        ax1.set_xticks(all_x, all_names, fontsize=10, rotation=90)

    # Horizontal lines (only if there is data range)
    if len(all_x) > 0:
        ax1.hlines(0, 0, len(all_names), lw=1, ls='--', color="k")
        ax1.hlines(dG_threshold, 0, len(all_names), lw=2, ls='--', color="b", alpha=0.5)

    ax1.set_ylabel(r"$\Delta G_{sol}$ [kJ/mol]", fontsize=14)
    ax1.legend(fontsize=12)
    ax1.set_title("Hydrophobic/aromatic sites (+ve solvation free energy)")

    # =====================================================
    # Panel 2: H-bond Acceptors (protein as acceptor)
    # =====================================================
    Ac1, Ac2, Ac_name1, Ac_name2 = [], [], [], []
    for key, values in data['hbonds']['acceptor'].items():
        resid = int(key)
        if resid in Index_frag1 and values['count'] > 5:
            Ac1.append(values['count']); Ac_name1.append(values['name'])
        elif resid in Index_frag2 and values['count'] > 5:
            Ac2.append(values['count']); Ac_name2.append(values['name'])

    # Safe x positions
    x1_ac = np.arange(len(Ac1))
    x2_ac = np.arange(len(Ac2)) + (x1_ac[-1] + 1 if len(x1_ac) > 0 else 0)

    if len(Ac1) > 0:
        ax2.bar(x1_ac, Ac1, label="Chain A", edgecolor="k", lw=1, width=0.5)
    if len(Ac2) > 0:
        ax2.bar(x2_ac, Ac2, label="Chain B", edgecolor="k", lw=1, width=0.5, fc='tab:red')

    all_ac_names = Ac_name1 + Ac_name2
    all_ac_x = np.concatenate([x1_ac, x2_ac]) if len(all_ac_names) > 0 else []
    if len(all_ac_x) > 0:
        ax2.set_xticks(all_ac_x, all_ac_names, rotation=90)
        ax2.hlines(donor_threshold, 0, len(all_ac_names), lw=2, ls='--', color="b", alpha=0.5)

    ax2.set_ylabel("H-bond/(frame*sasa)", fontsize=14)
    ax2.legend(fontsize=12)
    ax2.set_title("Donor pharmacophore site (i.e acceptors from the protein)")

    # =====================================================
    # Panel 3: H-bond Donors (protein as donor)
    # =====================================================
    Do1, Do2, Do_name1, Do_name2 = [], [], [], []
    for key, values in data['hbonds']['donor'].items():
        resid = int(key)
        if resid in Index_frag1 and values['count'] > 1:
            Do1.append(values['count']); Do_name1.append(values['name'])
        elif resid in Index_frag2 and values['count'] > 1:
            Do2.append(values['count']); Do_name2.append(values['name'])

    x1_do = np.arange(len(Do1))
    x2_do = np.arange(len(Do2)) + (x1_do[-1] + 1 if len(x1_do) > 0 else 0)

    if len(Do1) > 0:
        ax3.bar(x1_do, Do1, label="Chain A", edgecolor="k", lw=1, width=0.3)
    if len(Do2) > 0:
        ax3.bar(x2_do, Do2, label="Chain B", edgecolor="k", lw=1, width=0.4, fc='tab:red')

    all_do_names = Do_name1 + Do_name2
    all_do_x = np.concatenate([x1_do, x2_do]) if len(all_do_names) > 0 else []
    if len(all_do_x) > 0:
        ax3.set_xticks(all_do_x, all_do_names, rotation=90)
        ax3.hlines(acceptor_threshold, 0, len(all_do_names), lw=2, ls='--', color="b", alpha=0.5)

    ax3.set_ylabel("H-bond/(frame*sasa)", fontsize=14)
    ax3.legend(fontsize=12)
    ax3.set_title("Acceptors pharmacophore site (i.e donors from the protein)")

    # =====================================================
    # Panel 4: Ion Sites
    # =====================================================
    Cat1, Cat2, Cat_name1, Cat_name2 = [], [], [], []
    An1, An2, An_name1, An_name2 = [], [], [], []
    
    for (key_cat, val_cat), (key_an, val_an) in zip(
            data['ions']['cation'].items(),
            data['ions']['anion'].items()):
    
        if val_cat['count'] == 0.0 and val_an['count'] == 0:
            continue
    
        if int(key_cat) in Index_frag1:
            Cat1.append(val_cat['count'])
            Cat_name1.append(val_cat['name'])
            An1.append(val_an['count'])
            An_name1.append(val_an['name'])
        else:
            Cat2.append(val_cat['count'])
            Cat_name2.append(val_cat['name'])
            An2.append(val_an['count'])
            An_name2.append(val_an['name'])
    
    # X positions
    x1_cat = np.arange(len(Cat1))
    x2_cat = np.arange(len(Cat2)) + (x1_cat[-1] + 1 if len(x1_cat) > 0 else 0)
    
    # Bars for chain A
    if len(Cat1) > 0:
        ax4.bar(x1_cat, Cat1, width=0.3, label="Na$^+$", edgecolor="k", fc="tab:blue")
        ax4.bar(x1_cat + 0.33, An1, width=0.3, label="Cl$^-$", edgecolor="k", fc="tab:orange")
    
    # Bars for chain B (handle label duplication)
    if len(Cat2) > 0:
        na_label = None if len(Cat1) > 0 else "Na$^+$"
        cl_label = None if len(Cat1) > 0 else "Cl$^-$"
        ax4.bar(x2_cat, Cat2, width=0.3, label=na_label, edgecolor="k", fc="tab:blue")
        ax4.bar(x2_cat + 0.33, An2, width=0.3, label=cl_label, edgecolor="k", fc="tab:orange")
    
    # X ticks
    all_ion_names = An_name1 + An_name2
    all_ion_x = np.concatenate([x1_cat, x2_cat]) if len(all_ion_names) > 0 else []
    if len(all_ion_x) > 0:
        ax4.set_xticks(all_ion_x + 0.15, all_ion_names, rotation=90)
    
    # Threshold highlighting (difference >= ion_threshold)
    if len(Cat1) > 0 and len(An1) > 0:
        diff1 = np.abs(np.array(Cat1) - np.array(An1))
        index_1 = np.where(diff1 >= ion_threshold)[0]
        for i in index_1:
            ax4.bar(x1_cat[i] + 0.12, max(An1[i], Cat1[i]), width=0.6,
                    edgecolor="g", lw=2, fc=(1, 1, 1, 0))
    
    if len(Cat2) > 0 and len(An2) > 0:
        diff2 = np.abs(np.array(Cat2) - np.array(An2))
        index_2 = np.where(diff2 >= ion_threshold)[0]
        for i in index_2:
            ax4.bar(x2_cat[i] + 0.12, max(An2[i], Cat2[i]), width=0.6,
                    edgecolor="g", lw=2, fc=(1, 1, 1, 0))
    
    # Background spans (colored regions for chains)
    if len(An1) > 0:
        ax4.axvspan(-0.1, len(An1) - 0.3, color='blue', alpha=0.3, zorder=0)
    if len(An2) > 0:
        start = len(An1) - 0.3 if len(An1) > 0 else -0.1
        ax4.axvspan(start, start + len(An2), color='red', alpha=0.3, zorder=0)
    
    # Set xlim to avoid cut-off bars
    total_bars = len(Cat1) + len(Cat2)
    if total_bars > 0:
        ax4.set_xlim(-0.5, total_bars - 0.2)
    
    # Labels & annotations
    ax1.set_ylabel(r"$\mathrm{\Delta G_{sol} ~[kJ/mol]}$", fontsize=12)
    ax2.set_ylabel("H-bond donor", fontsize=12)
    ax3.set_ylabel("H-bond acceptor", fontsize=12)
    ax4.set_ylabel("Ion bound frames (%)", fontsize=12)
    
    ax1.annotate("(A)", xy=(0., 1.03), xycoords="axes fraction", fontsize=14)
    ax2.annotate("(B)", xy=(0., 1.03), xycoords="axes fraction", fontsize=14)
    ax3.annotate("(C)", xy=(0., 1.03), xycoords="axes fraction", fontsize=14)
    ax4.annotate("(D)", xy=(0., 1.03), xycoords="axes fraction", fontsize=14)
    
    if len(Cat1) > 0:
        ax4.annotate("Chain A", xy=(0.08, 0.6), xycoords="axes fraction", fontsize=12, color="tab:blue")
    if len(Cat2) > 0:
        ax4.annotate("Chain B", xy=(0.7, 0.6), xycoords="axes fraction", fontsize=12, color="tab:red")
    
    ax4.legend()
    
    # ---------------- Save or Show ---------------- #
    if out_file:
        plt.savefig(out_file, dpi=300, bbox_inches="tight")
    else:
        plt.show()


if __name__ == "__main__":
    main()
