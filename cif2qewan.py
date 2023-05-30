#!/usr/bin/env python
"""
Usage:
  cif2qewan.py <cif_file> <toml_file> [--so] [--mag]

Options:
  --so     input including spin-orbit couplings
  --mag    ferromagnetic calculations
"""

from docopt import docopt
import sys
import os
import re
import itertools
import numpy as np
import toml

class qe_wannier_in:
    def __init__(self, cif_file, toml_file, so, mag):

        info = toml.load(toml_file)
        self.cif2cell_path=info["cif2cell_path"]
        self.scf_k_resolution = info["scf_k_resolution"]
        self.degauss = info["degauss"]
        self.pseudo_dir = info["pseudo_dir"]
        self.pp_file_name = info["pp_file_name"]
        self.info_pw2wan = info["pw2wan"]

        self.so = so
        self.mag = mag

        self.lines = self.cif2cell_scf_in(cif_file)

        self.system_str = "&system\n"
        self.control_str = "&control\n"
        self.electrons_str = "&electrons\n"

        ntyp, nat = self.read_set_system()
        ecut_rho, ecut_wfc = self.read_set_pseudo_other(ntyp, nat, self.pp_file_name)

        self.set_system2(ntyp, ecut_rho, ecut_wfc)

        self.set_control()

        self.set_electrons(conv_thr="1.0d-8")


    def cif2cell_scf_in(self, cif_file):
        # call cif2cell and return result file (lines)
        cif_scf_in = "cif_scf.in"
        if(not os.path.exists(cif_scf_in)):
            # replace "(" and ")"
            cif_file = sys.argv[1].replace("(", "\(").replace(")", "\)")
            # call cif2cell
            os.system(self.cif2cell_path + " -p pwscf --setup-all --k-resolution={:0.3f} --print-digits=10 -o {} {}".format(self.scf_k_resolution, cif_scf_in, cif_file))
        return open(cif_scf_in).readlines()

    def read_set_system(self):
        # set ibrav, A, ntyp, nat
        # return ntyp, nat
        for line in self.lines:
            if ("ibrav" in line): self.system_str += line
            if ("A =" in line): 
                self.system_str += line
                self.alat = float(line.split("=")[1])
            if ("ntyp =" in line): 
                self.system_str += line
                ntyp = int(line.split("=")[1])
            if ("nat =" in line):
                self.system_str += line
                nat = int(line.split("=")[1])
        return ntyp, nat

    def set_system2(self, ntyp, ecut_rho, ecut_wfc):
        self.system_str += "  ecutwfc = {}\n".format(ecut_wfc)
        self.system_str += "  ecutrho = {}\n".format(ecut_rho)
        self.system_str += "  occupations = 'smearing'\n"
        self.system_str += "  smearing = 'm-p'\n"
        self.system_str += "  degauss = {:0.3f}\n".format(self.degauss)
        if(self.mag):
            self.system_str += "  nspin = 2\n"
            for i in range(ntyp):
                self.system_str += "  starting_magnetization(" + str(i+1) + ") = 3.0\n"
        elif(self.so):
            self.system_str += "  lspinorb = .true.\n"
            self.system_str += "  noncolin = .true.\n"

    def set_control(self):
        self.control_str += "  calculation = 'scf'\n"
        self.control_str += "  restart_mode = 'from_scratch'\n"
        self.control_str += "  prefix = 'pwscf'\n"
        self.control_str += "  tstress = .true.\n"
        self.control_str += "  tprnfor = .true.\n"
        self.control_str += "  pseudo_dir = '{}'\n".format(self.pseudo_dir)
        self.control_str += "  outdir = './work'\n"
        self.control_str += "  wf_collect = .true.\n"
        self.control_str += "  disk_io = 'low'\n"

    def set_electrons(self, conv_thr):
        self.electrons_str += "  mixing_mode = 'plain'\n"
        self.electrons_str += "  mixing_beta = 0.1\n"
        self.electrons_str += "  conv_thr = {}\n".format(conv_thr)

    def read_set_pseudo_other(self, ntyp, nat, pp_file_name):
        ecut_rho = 0
        ecut_wfc = 0
        pslist = pseudo_list(self.pseudo_dir, self.pp_file_name)
        self.pseudo_str = "ATOMIC_SPECIES\n"
        self.projection_str = ""
        num_wann_list = {}
        nexclude_list = {}
        self.atom_list = []
        self.atom_pos_list = []
        self.num_wann = 0
        self.nexclude = 0
        for i, line in enumerate(self.lines):
            if("CELL_PARAM" in line):
                self.cellparam_str = "".join(self.lines[i:i+4])
                self.a1 = np.array([ float(x) for x in self.lines[i+1].split() ])
                self.a2 = np.array([ float(x) for x in self.lines[i+2].split() ])
                self.a3 = np.array([ float(x) for x in self.lines[i+3].split() ])
            if ("ATOMIC_POSITIONS" in line):
                self.atompos_str = "".join(self.lines[i:i+nat+1])
                self.wan_atompos_str = "".join(self.lines[i+1:i+nat+1])
                for j in range(nat):
                    atm = self.lines[i+j+1].split()[0]
                    self.atom_list.append(atm)
                    self.atom_pos_list.append([float(x) for x in (self.lines[i+j+1].split()[1:4])])
                    self.num_wann += num_wann_list[atm]
                    self.nexclude += nexclude_list[atm]
            if ("K_POINTS" in line): 
                self.kpoints_str = "".join(self.lines[i:i+2])
                self.kmesh = [ int(x) for x in self.lines[i+1].split()[0:3] ]
            if("ATOMIC_SPECIES" in line):
                for j in range(ntyp):
                    line = self.lines[i+j+1]
                    atm = line[:5].strip()
                    ps = pslist.pseudo(atm)
                    ecut_rho = max(ecut_rho, ps.ecut_rho())
                    ecut_wfc = max(ecut_wfc, ps.ecut_wfc())
                    self.pseudo_str += re.sub('[A-Za-z]+_PSEUDO', ps.pseudo_file(), line)
                    if(self.so and not self.mag):
                        self.pseudo_str = self.pseudo_str.replace('.pbe', '.rel-pbe')
                    if(ps.wannier_orb != ""):
                        self.projection_str += "{}:{}\n".format(atm, ",".join(list(ps.wannier_orb)))
                    num_wann_list[atm] = ps.num_wann
                    nexclude_list[atm] = ps.nexclude
        return ecut_rho, ecut_wfc

    def convert2nscf(self):
        self.control_str = self.control_str.replace("'scf'", "'nscf'")

        system_add_str  = "  nosym = .true.\n"
        if(self.so or self.mag):
            system_add_str += "  nbnd = {}\n".format((self.nexclude + self.num_wann*3)*2)
        else:
            system_add_str += "  nbnd = {}\n".format(self.nexclude + self.num_wann*3)
        self.system_str = self.system_str.replace("&system\n", "&system\n" + system_add_str)

        if(self.mag):
            if(self.so):
                mag_str = "  lspinorb = .true.\n"
            else:
                mag_str = "  lspinorb = .false.\n"
            mag_str += "  noncolin = .true.\n"
            mag_str += "  lforcet = .true.\n"
            mag_str += "  angle1 = 0\n"
            mag_str += "  angle2 = 0\n"
            self.system_str = self.system_str.replace("  nspin = 2\n", mag_str)

            if(self.so):
                self.pseudo_str = self.pseudo_str.replace('.pbe', '.rel-pbe')

        self.electrons_str = re.sub("  conv_thr.*\n", "  conv_thr = 1.d-10\n", self.electrons_str)
        #self.electrons_str += "  diago_full_acc = .true.\n"

        self.nscfk = [ min( max(nk, 4), 8 ) for nk in self.kmesh ]
        self.kpoints_str = "K_POINTS {crystal}\n"
        self.kpoints_str += "{}\n".format(np.prod( self.nscfk ))
        self.wan_kmesh = ""
        for kx, ky, kz in itertools.product( range(self.nscfk[0]), range(self.nscfk[1]), range(self.nscfk[2]) ):
            kstr = "{:15.10f} {:15.10f} {:15.10f} {:15.10f}\n".format( float(kx)/self.nscfk[0], float(ky)/self.nscfk[1], float(kz)/self.nscfk[2], 1.0/np.prod(self.nscfk) )
            self.kpoints_str += kstr
            self.wan_kmesh += kstr

    def shift_k_nscf(self):
        """
        You MUST call this function after convert2nscf() was called.
        Creating nscf to calculate at shifted k-points from the original k-points.
        The shifted value is a half of one k-mesh in each of x,y,z directions.
        """
        self.control_str += "  verbosity = 'high'\n"
        self.system_str = self.system_str.replace("  nosym = .true.\n", "")
        if(self.so or self.mag):
            nbnd = (self.nexclude + int(self.num_wann*1.5))*2
        else:
            nbnd = self.nexclude + int(self.num_wann*1.5)
        self.system_str = re.sub("  nbnd.*\n", "  nbnd = {}\n".format(nbnd), self.system_str)
        self.electrons_str = self.electrons_str.replace("  diago_full_acc = .true.\n", "")
        self.electrons_str = re.sub("  conv_thr.*\n", "  conv_thr = 1.d-8\n", self.electrons_str)
        self.kpoints_str = "K_POINTS {automatic}\n"
        self.kpoints_str += "{0[0]} {0[1]} {0[2]}  1 1 1\n".format(self.nscfk)

    def convert2band(self):
        self.control_str = self.control_str.replace("'nscf'", "'bands'")
        self.kpoints_str = "K_POINTS {crystal_b}\n"
        self.kpoints_str += "{}\n".format(len(self.tick_labels) - self.tick_labels.count(""))
        for i in range(len(self.tick_labels)):
            if(i != (len(self.tick_labels)-1) and self.tick_labels[i+1] == ""):
                kstr = "{:15.10f} {:15.10f} {:15.10f}     {}    !  {}\n".format(self.tick_locs[i][0], self.tick_locs[i][1], self.tick_locs[i][2], 0, self.tick_labels[i])
                self.kpoints_str += kstr
            elif(self.tick_labels[i] != ""):
                kstr = "{:15.10f} {:15.10f} {:15.10f}    {}    !  {}\n".format(self.tick_locs[i][0], self.tick_locs[i][1], self.tick_locs[i][2], 20, self.tick_labels[i])
                self.kpoints_str += kstr

    def calc_bands_seekpath(self):
        try:
            import seekpath

        except ImportError:
            print("Failed to import seek path. Simple kpath is used instead.")
            self.tick_labels = ['R', 'G', 'X', 'M', 'G']
            self.tick_locs = [[0.5, 0.5, 0.5], [0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [0.5, 0.5, 0.0], [0.0, 0.0, 0.0]]
            return

        import pymatgen as mg

        cell = np.array([self.a1, self.a2, self.a3])
        pos = self.atom_pos_list
        z = [mg.Element(s).Z for s in self.atom_list]
        kpath = seekpath.getpaths.get_explicit_k_path([cell, pos, z])

        new_b = kpath['reciprocal_primitive_lattice']
        m = np.matmul(new_b,cell.T) / (2 * np.pi)
        self.kpoints_rel = [ np.matmul(k, m) for k in kpath["explicit_kpoints_rel"] ]

        kpoints_labels = kpath["explicit_kpoints_labels"]

        self.tick_locs = []
        self.tick_labels = []

        for i, label in enumerate(kpoints_labels):
            if(label == ""): continue
            label = label.replace("GAMMA","G")
            label = label.replace("SIGMA","S")
            if(i != 0 and kpoints_labels[i-1] != ""):
                self.tick_labels.extend(["", label])
                self.tick_locs.extend([np.array([0.0, 0.0, 0.0]), self.kpoints_rel[i]])
            else:
                self.tick_labels.append(label)
                self.tick_locs.append(self.kpoints_rel[i])

    def write_pwscf_in(self, pwscf_in):
        with open(pwscf_in, "w") as fp:
            fp.write(self.control_str + "/\n")
            fp.write(self.system_str + "/\n")
            fp.write(self.electrons_str + "/\n")
            fp.write(self.cellparam_str + "\n")
            fp.write(self.pseudo_str + "\n")
            fp.write(self.atompos_str + "\n")
            fp.write(self.kpoints_str + "\n")

    def write_band_in(self, band_in):
        with open(band_in, "w") as fp:
            fp.write("&bands\n")
            fp.write(" prefix = 'pwscf'\n")
            fp.write(" outdir = './work/'\n")
            fp.write(" filband = 'bands.out'\n")
            fp.write("/\n")
        
    def write_pw2wan(self, pw2wan):
        info = self.info_pw2wan
        with open(pw2wan, "w") as fp:
            fp.write("&inputpp\n")
            fp.write(" outdir = './work'\n")
            fp.write(" prefix = 'pwscf'\n")
            fp.write(" seedname = 'pwscf'\n")
            fp.write(" spin_component = 'none'\n")
            fp.write(" write_mmn = .true.\n")
            fp.write(" write_amn = .true.\n")
            fp.write(" write_unk = {}\n".format(info["write_unk"]))
            fp.write("/\n")

    def write_wannier(self, wannier_in):
        with open(wannier_in, "w") as fp:
            fp.write("! generated by {}\n".format(__file__.split("/")[-1]))
            so_factor = 1
            if(self.so or self.mag): so_factor = 2
            fp.write("num_bands = {}\n".format(self.num_wann*3*so_factor))
            fp.write("num_wann  = {}\n".format(self.num_wann*so_factor))
            if(self.nexclude > 0):
                fp.write("exclude_bands = 1-{}\n\n".format(self.nexclude*so_factor))

            fp.write("dis_num_iter = 200\n")
            fp.write("num_iter = 0\n\n")
            fp.write("dis_froz_max = -200\n")
            fp.write("dis_froz_min = -200\n\n")


            if(self.so or self.mag):
                fp.write("spinors = .true.\n\n")

            fp.write("begin projections\n")
            fp.write(self.projection_str)
            fp.write("end projections\n\n")

            fp.write("bands_plot = .true.\n")
            fp.write("write_hr = .true.\n")
            fp.write("write_tb = .true.\n")
            fp.write("fermi_surface_plot = .true.\n")
            fp.write("wannier_plot = .true.\n")
            fp.write("\n")

            fp.write("begin unit_cell_cart\n")
            fp.write("ang\n")
            fp.write("{0[0]:12.7f} {0[1]:12.7f} {0[2]:12.7f}\n".format(self.a1 * self.alat))
            fp.write("{0[0]:12.7f} {0[1]:12.7f} {0[2]:12.7f}\n".format(self.a2 * self.alat))
            fp.write("{0[0]:12.7f} {0[1]:12.7f} {0[2]:12.7f}\n".format(self.a3 * self.alat))
            fp.write("end unit_cell_cart\n\n")

            fp.write("begin atoms_frac\n")
            fp.write(self.wan_atompos_str)
            fp.write("end atoms_frac\n\n")

            fp.write("mp_grid: {0[0]} {0[1]} {0[2]}\n\n".format(self.nscfk))

            fp.write("begin kpoints\n")
            fp.write(self.wan_kmesh)
            fp.write("end kpoints\n\n")

            fp.write("Begin Kpoint_Path\n")
            for i in range(len(self.tick_labels) - 1):
                if(self.tick_labels[i] != "" and self.tick_labels[i+1] != ""):
                    fp.write("{0} {1[0]:14.10f} {1[1]:14.10f} {1[2]:14.10f}  {2} {3[0]:14.10f} {3[1]:14.10f} {3[2]:14.10f}\n".format(self.tick_labels[i], self.tick_locs[i], self.tick_labels[i+1], self.tick_locs[i+1]))
            fp.write("End Kpoint_Path\n")

    def write_proj(self, proj_in):
        with open(proj_in, "w") as fp:
            fp.write("&projwfc\n")
            fp.write(" prefix = 'pwscf'\n")
            fp.write(" outdir = './work'\n")
            fp.write(" kresolveddos = .false.\n")
            fp.write(" degauss = {:0.3f}\n".format(self.degauss))
            fp.write(" Emax = \n")
            fp.write(" Emin = \n")
            fp.write("/\n")


class pseudo_wannier:
    """
    pseudo potential including information for wannier
    """
    def __init__(self, pseudo_dir, atm, type, nexclude, wannier_orb):
        self.atm = atm
        self.type = type
        self.nexclude = nexclude
        self.wannier_orb = wannier_orb
        self._ecut_wfc = 0
        self._ecut_rho = 0
        self.pseudo_dir = pseudo_dir

        self.num_wann = 0
        for s in wannier_orb:
            if(s == "s"): self.num_wann += 1
            if(s == "p"): self.num_wann += 3
            if(s == "d"): self.num_wann += 5
            if(s == "f"): self.num_wann += 7

    def set_ecut(self):
        if(self._ecut_rho == 0 or self._ecut_wfc == 0):
            with open(self.pseudo_dir + "/" + self.pseudo_file(), "r") as fp:
                for line in fp.readlines():
                    if("Suggested minimum cutoff for wavefunctions" in line):
                        self._ecut_wfc = float(line.split()[5])
                    if("Suggested minimum cutoff for charge density" in line):
                        self._ecut_rho = float(line.split()[6])
                        break

    def ecut_wfc(self):
        self.set_ecut()
        return self._ecut_wfc

    def ecut_rho(self):
        self.set_ecut()
        return self._ecut_rho

    def pseudo_file(self):
        return "{}.{}.UPF".format(self.atm, self.type)

class pseudo_list:
    def __init__(self, pseudo_dir, pp_file_name):
        self.ps_list = {}
        pp_info = self.read_pp_info(pp_file_name)
        for key, values in pp_info.items():
            if values is not None:
                self.ps_list[key] = pseudo_wannier(pseudo_dir, key, values[0], values[1], values[2])
            else:
                self.ps_list[key] = None

    def pseudo(self, atm):
        return self.ps_list[atm]

    def read_pp_info(self, pp_file_name):
        import pandas as pd          
        csv_input = pd.read_csv(pp_file_name, sep=",")
        #print(csv_input)
        dict_pp = {}
        for value in csv_input.values:
            if value[1] is np.nan:
                dict_pp[value[0]] = None
            else:
                if value[2] is np.nan:
                    value[2] = 0
                elif value[3] is np.nan:
                    value[3] = ""
                dict_pp[value[0]]=(value[1], int(value[2]), value[3])
        return dict_pp

if __name__ == '__main__':
    args = docopt(__doc__)

    cif_file = args["<cif_file>"]
    toml_file = args["<toml_file>"]

    qe_wan = qe_wannier_in(cif_file, toml_file, args["--so"], args["--mag"])

    qe_wan.write_pwscf_in("scf.in")

    qe_wan.convert2nscf()

    qe_wan.write_pwscf_in("nscf.in")

    qe_wan.calc_bands_seekpath()

    qe_wan.write_pw2wan("pw2wan.in")

    qe_wan.write_wannier("pwscf.win")

    qe_wan.write_proj("proj.in")

    if not os.path.exists("check_wannier"): os.mkdir("check_wannier")
    qe_wan.shift_k_nscf()
    qe_wan.write_pwscf_in("check_wannier/nscf.in")

    if not os.path.exists("band"): os.mkdir("band")
    qe_wan.convert2band()
    qe_wan.write_pwscf_in("band/nscf.in")
    qe_wan.write_band_in("band/band.in")
