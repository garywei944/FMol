import networkx as nx
from rdkit.Chem import Descriptors, rdmolops, MolFromSmiles, AddHs, AllChem, MolToInchi
from abc import ABC, abstractmethod

from mcts import parameters as p
from mcts.properties import sascorer
from mcts.properties.dft import calcul_dft


def decorable(cls):
    """
    function to decorate a class with "<"

    :param cls: class to decorate
    :type cls: class
    :return: class
    """
    cls.__lshift__ = lambda property_decorated, function: function(property_decorated)
    return cls


@decorable
class Property(ABC):
    """
    Mother class of properties
    """
    def __init__(self, smiles):
        """
        Create the new property for the SMILES given in parameter

        :param smiles: SMILES to analyse
        :type smiles: SMILES
        """
        self.smiles = smiles

    @abstractmethod
    def calculate(self):
        """
        Abstract method to calculate the properties of a SMILES

        :return: None
        """
        pass


class Property2D(Property):
    """
    2D properties, based on RDKit representation
    """
    def __init__(self, smiles):
        super().__init__(smiles)

    def calculate(self):
        """
        Check if the SMILES is valid then update the info.

        :return: RDKit Mol object
        """
        try:
            m = MolFromSmiles("".join(p.config['long_prefix']) + "".join(self.smiles.element))
            self.smiles.properties[p.s_valid] = False
            if m is not None:
                m = AddHs(m)
                AllChem.EmbedMolecule(m)
                AllChem.UFFOptimizeMolecule(m)
                self.smiles.properties["InChI"] = MolToInchi(m)
        except Exception as e:
            print("Error rdkit : " + repr(e))
            m = None
        if m is not None:
            self.smiles.properties[p.s_valid] = True
            with p.lock_update_data:
                p.tree_info[p.info_good] += 1
                self.smiles.properties[p.s_id] = p.tree_info[p.info_good]
        else:
            with p.lock_update_data:
                p.tree_info[p.info_bad] += 1
                self.smiles.properties[p.s_id] = p.tree_info[p.info_bad]
        return m


class PropertyDecorator(Property):
    """
    Abstract class for property decorator
    """
    def __init__(self, property_decorated):
        """
        Create the new property decorator

        :param property_decorated: the property to decorate
        """
        super().__init__(property_decorated.smiles)
        self.property = property_decorated

    @abstractmethod
    def calculate(self):
        pass


class SAScoreProperty2DDecorator(PropertyDecorator):
    """
    Decorator for SA Score
    """
    def __init__(self, property_decorated):
        super().__init__(property_decorated)

    def calculate(self):
        """
        Call the calculate method of the property then calculate the SA score before returning the RDKit mol

        :return: RDKit mol
        """
        m = self.property.calculate()
        if m is not None:
            try:
                with p.lock_sa_score:
                    self.smiles.properties[p.s_sa] = sascorer.calculate_score(m)
            except Exception as e:
                print("Error SA : " + repr(e))
                print(e)
                self.smiles.properties[p.s_valid] = False
                m = None
        return m


class CycleProperty2DDecorator(PropertyDecorator):
    """
    Cycle property
    """
    def __init__(self, property_decorated):
        super().__init__(property_decorated)

    def calculate(self):
        """
        Call the calculate method of the property then calculate the size of the longest cycle
        before returning the RDKit mol

        :return: RDKit mol
        """
        m = self.property.calculate()
        if m is not None:
            try:
                cycle_list = nx.cycle_basis(nx.Graph(rdmolops.GetAdjacencyMatrix(m)))
                self.smiles.properties[p.s_cycle] = max([len(j) for j in cycle_list]) if cycle_list else 0
            except Exception as e:
                print("Error Cycle : " + repr(e))
                self.smiles.properties[p.s_valid] = False
                m = None
        return m


class LogPProperty2DDecorator(PropertyDecorator):
    """
    LogP property
    """
    def __init__(self, property_decorated):
        super().__init__(property_decorated)

    def calculate(self):
        """
        Call the calculate method of the property then calculate the LogP before returning the RDKit mol

        :return: RDKit mol
        """
        m = self.property.calculate()
        if m is not None:
            try:
                self.smiles.properties[p.s_logp] = Descriptors.MolLogP(m)
            except Exception as e:
                print("Error LogP : " + repr(e))
                self.smiles.properties[p.s_valid] = False
                m = None
        return m


class DFTPropertyDecorator(PropertyDecorator):
    """
    DFT property
    """
    def __init__(self, property_decorated):
        super().__init__(property_decorated)

    def calculate(self):
        """
        Call the calculate method of the property then launch the DFT before returning the RDKit mol

        :return: RDKit mol
        """
        m = self.property.calculate()
        if self.smiles.properties[p.s_valid]:
            self.smiles.properties[p.s_dft] = calcul_dft(self.smiles.properties[p.s_id],
                                                  "".join(p.config['long_prefix']) + "".join(self.smiles.element[:-1]),
                                                  m)
        return m


if __name__ == '__main__':
    from mcts.smiles import SMILES
    import threading
    p.lock_sa_score = threading.Lock()
    p.config["long_prefix"] = []
    s = SMILES(["Cl", "O", "C"])
    print(s.properties)
    pa = Property2D(s) << CycleProperty2DDecorator << LogPProperty2DDecorator << SAScoreProperty2DDecorator
    pa.calculate()
    print(s.properties)
