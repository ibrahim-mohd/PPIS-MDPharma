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
    
    ##################### Thresholds
    
    parser.add_argument('-acceptor_th', dest='acceptor_threshold', type=float, default=40, help='acceptor threshold')
    parser.add_argument('-donor_th', dest='donor_threshold', type=float, default=40, help='acceptor threshold')
    parser.add_argument('-dG_th', dest='dG_threshold', type=float, default=0.4, help='ΔG threshold,')
    parser.add_argument('-ion_th', dest='ion_threshold', type=float, default=30, help='difference between number of cations and anions')

    # Figure options
    parser.add_argument('-o',dest='out_file', type=str, default="features.png", help="Output image file (e.g. plot.png). If None, shows figure.")
    parser.add_argument('-figsize',dest='figsize', type=float, nargs=2, default=(12, 12), help="Figure size (width height)")
    
    return parser.parse_args()


def main():
    
    args               = parse_args()
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
    with open(pkl_file, "rb") as f: data = pickle.load(f)

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
    x1 = np.arange(len(dG1))
    x2 = np.arange(x1[-1] + 1, len(dG2) + len(x1))

    ax1.bar(x1, dG1, yerr=std1, label="Chain A", edgecolor="k", lw=1, capsize=5)
    ax1.bar(x2, dG2, yerr=std2, label="Chain B", edgecolor="k", lw=1, capsize=5, fc='tab:red')

    ax1.set_xticks(np.concatenate([x1, x2]), resname1 + resname2, fontsize=10, rotation=90)
    ax1.hlines(0, 0, len(resname1 + resname2), lw=1, ls='--', color="k")
    ax1.hlines(dG_threshold, 0, len(resname1 + resname2), lw=2, ls='--', color="b", alpha=0.5)

    ax1.set_ylabel(r"$\Delta G_{sol}$ [kJ/mol]", fontsize=14)
    ax1.legend(fontsize=12)
    ax1.set_title ("Hydrophobic/aromatic sites (+ve solvation free energy)")
    # =====================================================
    # Panel 2: H-bond Acceptors
    # =====================================================
    Ac1, Ac2, Ac_name1, Ac_name2 = [], [], [], []
    for key, values in data['hbonds']['acceptor'].items():
        resid = int(key)
        if resid in Index_frag1 and values['count'] > 5:
            Ac1.append(values['count']); Ac_name1.append(values['name'])
        elif resid in Index_frag2 and values['count'] > 5:
            Ac2.append(values['count']); Ac_name2.append(values['name'])

    x1 = np.arange(len(Ac1))
    x2 = np.arange(x1[-1] + 1, len(Ac2) + len(x1))

    ax2.bar(x1, Ac1, label="Chain A", edgecolor="k", lw=1, width=0.5)
    ax2.bar(x2, Ac2, label="Chain B", edgecolor="k", lw=1, width=0.5, fc='tab:red')

    ax2.set_xticks(np.concatenate([x1, x2]), Ac_name1 + Ac_name2, rotation=90)
    ax2.set_ylabel("H-bond/(frame*sasa)", fontsize=14)
    ax2.hlines(donor_threshold, 0, len(Ac_name1 + Ac_name2), lw=2, ls='--', color="b", alpha=0.5)
    ax2.legend(fontsize=12)
    ax2.set_title ("Donor pharmacophore site (i.e acceptors from the protein)")
    # =====================================================
    # Panel 3: H-bond Donors
    # =====================================================
    Do1, Do2, Do_name1, Do_name2 = [], [], [], []
    for key, values in data['hbonds']['donor'].items():
        resid = int(key)
        if resid in Index_frag1 and values['count'] > 1:
            Do1.append(values['count']); Do_name1.append(values['name'])
        elif resid in Index_frag2 and values['count'] > 1:
            Do2.append(values['count']); Do_name2.append(values['name'])

    x1 = np.arange(len(Do1))
    x2 = np.arange(x1[-1] + 1, len(Do2) + len(x1))

    ax3.bar(x1, Do1, label="Chain A", edgecolor="k", lw=1, width=0.3)
    ax3.bar(x2, Do2, label="Chain B", edgecolor="k", lw=1, width=0.4, fc='tab:red')

    ax3.set_xticks(np.concatenate([x1, x2]), Do_name1 + Do_name2, rotation=90)
    ax3.set_ylabel("H-bond/(frame*sasa)", fontsize=14)
    ax3.hlines(acceptor_threshold, 0, len(Do_name1 + Do_name2), lw=2, ls='--', color="b", alpha=0.5)
    ax3.legend(fontsize=12)
    ax3.set_title ("Acceptors pharmacophore site (i.e donors from the protein)")

    # =====================================================
    # Panel 4: Ion Sites
    # =====================================================
    Cat1, Cat2, Cat_name1, Cat_name2 = [], [], [], []
    An1, An2, An_name1, An_name2 = [], [], [], []

    for (key_cat, val_cat), (key_an, val_an) in zip(data['ions']['cation'].items(),
                                                    data['ions']['anion'].items()):
        if val_cat['count'] < 0.1 and val_an['count'] < 0.1:
            continue

        if int(key_cat) in Index_frag1:
            Cat1.append(val_cat['count']); Cat_name1.append(val_cat['name'])
            An1.append(val_an['count']);   An_name1.append(val_an['name'])
        else:
            Cat2.append(val_cat['count']); Cat_name2.append(val_cat['name'])
            An2.append(val_an['count']);   An_name2.append(val_an['name'])

    x1 = np.arange(len(Cat1))
    x2 = np.arange(x1[-1] + 1, len(Cat2) + len(x1))

    ax4.bar(x1, Cat1, width=0.3, label="Na$^+$, Chain A", edgecolor="k", fc="tab:blue")
    ax4.bar(x1 + 0.33, An1, width=0.3, label="Cl$^-$, Chain A", edgecolor="k", fc=(0, 0, 1, 0.5))
    ax4.bar(x2, Cat2, width=0.3, label="Na$^+$, Chain B", edgecolor="k", fc="tab:red")
    ax4.bar(x2 + 0.33, An2, width=0.3, label="Cl$^-$, Chain B", edgecolor="k", fc=(1, 0, 0, 0.5))

    ax4.set_xticks(np.concatenate([x1, x2]) + 0.15, An_name1 + An_name2, rotation=90)
    ax4.set_ylabel("Ion bound frame %age", fontsize=16)
    ax4.set_title ("Positive/negative ion binding sites")
    # -------------------- Threshold -------------------------
    index_1 = np.where ( np.abs (np.array (Cat1) - np.array (An1)) >= ion_threshold) [0]
    index_2 = np.where ( np.abs (np.array (Cat2) - np.array (An2)) >= ion_threshold) [0]

    if np.size (index_1):
        for i in index_1:
            ax4.bar(x1[i]+0.12, Cat1[i], width=0.6, edgecolor="g", lw=2, fc=(1,1,1,0))
           # ax4.bar(x1[i], An1[i]+0.33, width=0.3, edgecolor="g", lw=2, fc=(1,1,1,0))

    if np.size (index_2):
        for i in index_2:
            ax4.bar(x2[i]+0.12, Cat2[i], width=0.6, edgecolor="g", lw=2, fc=(1,1,1,0))
           # ax4.bar(x2[i], An2[i]+0.33, width=0.3, edgecolor="g", lw=2, fc=(1,1,1,0))
                
    # ---------------- Save or Show ---------------- #
    if out_file:
        plt.savefig(out_file, dpi=300, bbox_inches="tight")
    else:
        plt.show()


if __name__ == "__main__":
    main()
