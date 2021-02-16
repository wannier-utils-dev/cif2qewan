# coding : utf-8
# You can compare Wannier bands with PWSCF bands.
# Please use the same k-path and structure for calculation.

import numpy as np
import os.path
import matplotlib
matplotlib.use('Agg')

from matplotlib import pyplot as plt

scfout = "scf.out"
wannier_band = "pwscf_band.dat"
wannier_band_gnu = "pwscf_band.gnu"
pwscf_band = "./band/bands.out.gnu"
pwscf_win = "pwscf.win"


def get_ef_from_scfout():
    """
    Get Fermi Energy from QE scf.out.
    :param path_scf_out:
    :return: float(ef)
    """
    assert os.path.exists(scfout)

    ef = 0
    with open(scfout) as fp:
        for line in fp.readlines():
            if "Fermi" in line:
                ef = float(line.split()[-2])

    return ef


def get_band_data(band_file):
    """
    Get Band from qe(band.gnu)
    x(k-path) is  normalized.
    :return: x[nk], y[nk, nband]
    """

    data = np.loadtxt(band_file)
    x = data[:, 0]
    nband = np.sum(x == 0)
    nk = int(len(x) / nband)
    x = x[0:nk]/max(x)
    y = np.transpose(data[:, 1].reshape(nband, nk))

    return x, y


def get_froz_max():
    dis_froz_max = -200.0
    with open(pwscf_win) as fp:
        for line in fp.readlines():
            if "dis_froz_max" in line:
                dis_froz_max = float((line.split())[-1])

    return dis_froz_max


def get_klabel():
    """
    Get xtics from wannier_band.gnu
    :return: x(kpath)[num_klabel], label[num_klabel]
    """

    xtics = ""
    with open(wannier_band_gnu) as fp:
        for line in fp.readlines():
            if "set xtics" in line:
                xtics = line.replace("set xtics", "").strip()
                break

    xtics = xtics[1:-1].split(",")
    x_list = []
    label_list = []
    for xtic in xtics:
        label, pos = xtic.split()
        x_list.append(float(pos))
        label_list.append(label[1:-1].replace("G", "$\\Gamma$").replace("S", "$\\Sigma$"))

    x_list = np.array(x_list) / max(x_list)
    return x_list, label_list


def main():
    x, y = get_band_data(wannier_band)
    x_qe, y_qe = get_band_data(pwscf_band)

    ef = get_ef_from_scfout()

    froz_max = get_froz_max()

    klabel = get_klabel()

    # plot with matplotlib.pyplot
    plt.rcParams["font.size"] = 16

    plt.title("Red=QE, Black=Wannier")
    plt.plot(x_qe, y_qe - ef, c="r", lw=1, label="QE")
    plt.plot(x, y - ef, c="k", lw=0.5, label="Wannier")
    plt.ylabel("$E - E_{\mathrm{F}}$[eV]")
    plt.xticks(klabel[0], klabel[1])
    plt.xlim([0, 1])
    y_min = np.min(y-ef)
    y_max = np.max(y-ef)
    plt.ylim([y_min - 0.05*(y_max-y_min), y_max + 0.05*(y_max-y_min)])

    plt.savefig("./band/band_compare.png", bbox_inches='tight')
    plt.savefig("./band/band_compare.eps", bbox_inches='tight')


if __name__ == "__main__":
    main()
