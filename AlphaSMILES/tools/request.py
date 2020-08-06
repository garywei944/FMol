import json
import os

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw

from mcts import parameters as p
from tools.plot_wavelength import plot_wl


def select(data, starting_with='', wl_min=0, wl_max=float('inf'), unit="nm", f_min=0.0):
    """
    Select SMILES in data with oscillator strength peaks between bounds

    :param data: dict to analyse
    :type data: dict(SMILES)
    :param starting_with: beginning of the SMILES (you have to add an extra ' before the SMILES)
    :type starting_with: str
    :param wl_min: wavelength minimum for the search
    :type wl_min: float
    :param wl_max: wavelength maximum for the search
    :type wl_max: float
    :param unit: unit of the wavelength ( "nm" , "ev" or "cm-1" )
    :type unit: str
    :param f_min: oscillator strength minimum for the search
    :type f_min: float
    :return: list of selected smiles
    """
    selected_smiles = set()
    for smiles in data.keys():
        if data[smiles]['valid']:
            if smiles.startswith(starting_with):
                for line in data[smiles][p.s_dft]:
                    if wl_min <= line[unit] <= wl_max and line['f'] >= f_min:
                        selected_smiles.add(smiles)
                        print(smiles)
                        print(line)
    return selected_smiles


def smiles_to_image(i, smiles):
    """
    Convert a SMILES into image,
    the image name is "id" + "_2D.png"

    :param i: id of the SMILES given in parameter
    :type i: int
    :param smiles: String corresponding to the SMILES to convert to image
    :type smiles: str
    :return:
    """
    m = Chem.MolFromSmiles(smiles)
    AllChem.Compute2DCoords(m)
    Draw.MolToFile(m, "../data_out/" + p.config["configuration_name"] + "/plot/" + str(i) + '_2D.png')


def find(config_name, starting_with='', wl_min=0, wl_max=float('inf'), unit="nm", f_min=0.0, plot_wavelength=True):
    """
    Look for the SMILES with the given parameters

    :param config_name: name of the configuration to use
    :type config_name: str
    :param starting_with: the begging of the SMILES
    :type starting_with: str
    :param wl_min: wavelength minimum for the search
    :type wl_min: float
    :param wl_max: wavelength maximum for the search
    :type wl_max: float
    :param unit: unit of the wavelength ( "nm" , "ev" or "cm-1" )
    :type unit: str
    :param f_min: oscillator strength minimum for the search
    :type f_min: float
    :param plot_wavelength: if you want to get a plot of the wavelength or not
    :type plot_wavelength: bool
    :return: None
    """
    with open('../data_out/' + config_name + '/data.json') as d:
        data = json.load(d)

    with open('../mcts/configurations/' + config_name + ".json") as c:
        p.config = json.load(c)

    selected = select(data, starting_with, wl_min, wl_max, unit, f_min)

    if not os.path.isdir("../data_out/" + p.config["configuration_name"] + "/plot/"):
        os.mkdir("../data_out/" + p.config["configuration_name"] + "/plot/")

    for s in selected:
        smi = s[1:-1] if s[0] == "'" else s
        smi = "".join(p.config['long_prefix']) + smi
        smiles_to_image(data[s]['id'], smi)
        if plot_wavelength:
            plot_wl(data, s)


if __name__ == '__main__':
    find('generated', starting_with="'c1", wl_min=500, unit="nm", f_min=0.4)

