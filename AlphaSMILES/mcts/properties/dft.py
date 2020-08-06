from __future__ import absolute_import, division, print_function

import lzma
import os
import shutil
import time
from builtins import (str, open, int)

import cclib
from rdkit import Chem

from mcts import parameters as p


def delete_file(file):
    """
    Delete the file if exist

    :param file: the file to delete
    :type  file: str
    :return: None
    """
    if os.path.isfile(file):
        os.remove(file)


def compress_file(file):
    """
    Compress the file if exist, delete it and move it

    :param file: the file to compress and delete
    :type  file: str
    :return: None
    """
    if os.path.isfile(file):
        with lzma.open(file + ".tar.xz", "w") as f:
            f.write(open(file, 'rb').read())
        delete_file(file)
        move_file(file + ".tar.xz")


def move_file(file):
    """
    Move the file to generated/dft/

    :param file: the file to move
    :type  file: str
    :return: None
    """
    if os.path.isfile(file):
        shutil.move(file, p.r_dft + file)


def calcul_dft(id_smiles, smiles, m):
    """
    Convert the SMILES into RDKit mol,
    get the position of the atoms with obabel,
    launch OTP with gaussian
    if the OTP has a normal termination then launch TD
    return the oscillator strength and wavelength

    :param id_smiles: id of the SMILES
    :type id_smiles: int
    :param smiles: SMILES to analyse
    :type smiles: SMILES
    :param m: mol version of the SMILES
    :type m: RDKit mol
    :return: Oscillator Strength and wavelength in nm, ev and cm-1
    """
    print("Starting DFT on " + smiles)

    n_mol = str(id_smiles) + ".mol"
    n_xyz = str(id_smiles) + "_xyz"
    n_opt = str(id_smiles) + "_OPT.inp"
    n_opt_log = str(id_smiles) + "_OPT.log"
    n_td = str(id_smiles) + "_TD.inp"
    n_td_log = str(id_smiles) + "_TD.log"

    dft = []
    try:
        with open(n_mol, "w") as mol:
            mol.write(Chem.MolToMolBlock(m))
    except Exception as e:
        print(e)
        delete_file(n_mol)
        return dft
    command_obabel = "obabel -imol " + n_mol + " -oxyz -O " + n_xyz
    os.system(command_obabel)
    with open(n_xyz, "r") as xyz:
        position = ""
        for i, l in enumerate(xyz):
            if i >= 2:
                position += l

    delete_file(n_mol)
    delete_file(n_xyz)

    """
    TODO
    Code to use for GAUSSIAN
    I don't understand why but I get an 
    TypeError: slice indices must be integers or None or have an __index__ method
    
    while trying to open the log file with cclib...
    data = cclib.io.ccread(n_log)
    """

#     n_log = str(id_smiles) + "_DFT.log"
#     inp_file = """%Chk={id_smiles}_DFT
# %NProcShared={nb_core_dft}
# %mem=1200MB
# #P B3LYP/3-21G* opt gfprint pop=(full,HirshfeldEE)
#
# {id_smiles}_DFT
#
# 0 1
# {position}
#
#
# --Link1--
# %Chk={id_smiles}_DFT
# %NProcShared={nb_core_dft}
# %mem=1200MB
# #p B3LYP chkbasis guess=read geom=allcheck TD=(nstates=20) Gfprint
#
#
#
# """.format(id_smiles=id_smiles, nb_core_dft=p.config["nb_core_dft"], position=position)
#
#     n_inp = str(id_smiles) + "_DFT.inp"
#
#     with open(n_inp, "w") as inp:
#         inp.write(inp_file)
#
#     command_opt = "./mcts/properties/dft.sh " + n_inp
#     start = time.time()
#     os.system(command_opt)
#     stop = time.time()
#     print("Execution time DFT: " + repr(int(stop - start)) + "s")
#
#     with open(n_log, "r") as log:
#         last_line = log.readlines()[-1]
#
#     if "Normal termination" in last_line:
#         with open(n_log, "r") as log:
#             data = cclib.io.ccread(n_log)
#             i = 0
#             for line in log:
#                 if "Excited State" in line:
#                     val = line.split()[4:-1]
#                     dft.append(dict({"ev": float(val[0]),
#                                      "nm": float(val[2]),
#                                      "cm-1": float(data.etenergies[i]),
#                                      "f": float(val[-1].split("=")[1])}))
#                     i += 1
#
#     delete_file(n_inp)

    # Create inp file for OPT
    with open(n_opt, "w") as inp:
        inp.write("%Chk=" + str(id_smiles) + "\n")
        inp.write("%NProcShared=" + str(p.config["nb_core_dft"]) + "\n")
        inp.write("%mem=1200MB\n")
        inp.write("#P B3LYP/3-21G* opt gfprint pop=(full,HirshfeldEE)\n")
        inp.write("\n" + str(id_smiles) + "\n\n")
        inp.write("0 1\n")
        inp.write(position + "\n\n\n")

    # Calculate OPT
    command_opt = "./mcts/properties/dft.sh " + n_opt
    print("Starting OPT")
    start = time.time()
    os.system(command_opt)
    stop = time.time()
    print("Execution time OPT: " + repr(int(stop - start)) + "s")

    delete_file(n_opt)

    with open(n_opt_log, "r") as log:
        last_line = log.readlines()[-1]

    # if the OTP end up well
    if "Normal termination" in last_line:
        # Create inp file for TD
        with open(n_td, "w") as inp:
            inp.write("%Chk=" + str(id_smiles) + "\n")
            inp.write("%NProcShared=" + str(p.config["nb_core_dft"]) + "\n")
            inp.write("%mem=1200MB\n")
            inp.write("#p B3LYP chkbasis guess=read geom=allcheck TD=(nstates=20) Gfprint\n\n\n")

        start = time.time()
        # Calculate TD
        command_td = "./mcts/properties/dft.sh " + n_td
        print("Starting TD")
        os.system(command_td)
        stop = time.time()
        print("Execution time TD: " + repr(int(stop - start)) + "s")

        delete_file(n_td)

        with open(n_td_log, "r") as log:
            last_line = log.readlines()[-1]

        # if the TD end up well
        if "Normal termination" in last_line:
            # create a list of dict with the result of DFT
            with open(str(id_smiles) + "_TD.log", "r") as log:
                data = cclib.io.ccread(n_td_log)
                i = 0
                for line in log:
                    if "Excited State" in line:
                        val = line.split()[4:-1]
                        dft.append(dict({"ev": float(val[0]),
                                         "nm": float(val[2]),
                                         "cm-1": float(data.etenergies[i]),
                                         "f": float(val[-1].split("=")[1])}))
                        i += 1
    compress_file(n_opt_log)
    compress_file(n_td_log)

    delete_file(str(id_smiles) + ".chk")

    return dft

