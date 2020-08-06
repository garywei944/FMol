import json
import threading
import os

from mcts import parameters as p


def load_data_base():
    """
    Prepare the data_base

    :return: None
    """
    p.lock_access_data_base = threading.Lock()
    data_file = "data_out/" + p.config["data_base"] + ".json"
    if os.path.isfile(data_file):
        print("Data base loaded")
        print(data_file)
        with open(data_file, 'r') as f:
            p.data_base = json.load(f)
    else:
        print("Creation of a new data base")
        p.data_base = dict()


def save_data_base():
    """
    Save the data base in the "p.config["data_base"].json" file

    :return: None
    """
    with open("data_out/" + p.config["data_base"] + ".json", 'w') as f:
        json.dump(p.data_base, f)


def create(smiles, properties):
    """
    Add a new SMILES in the data base if it doesn't already exist.

    :param smiles: the SMILES to create
    :type smiles: str
    :param properties: The properties of the SMILES
    :type properties: set
    :return: True if the creation end up well, False otherwise
    """
    with p.lock_access_data_base:
        if smiles not in p.data_base:
            p.data_base[smiles] = properties
            return True
        else:
            return False


def select(smiles):
    """
    Select a SMILES from the data base, return None if the SMILES is not in the data base.

    :param smiles: SMILES to look for in the data base
    :type smiles: str
    :return: None or SMILES
    """
    with p.lock_access_data_base:
        return p.data_base.get(smiles, None)


def update(smiles, properties):
    """
    Update a SMILES with properties given in parameters.

    :param smiles: the SMILES to update
    :type smiles: str
    :param properties: The new properties
    :type properties: set
    :return: True if the update end up well, False otherwise
    """
    with p.lock_access_data_base:
        if smiles in p.data_base:
            p.data_base[smiles] = properties
            return True
        else:
            return False


def delete(smiles):
    """
    Delete the SMILES from the data base

    :param smiles: the SMILES to delete
    :type smiles: str
    :return: the properties of the SMILES if the deletion end up well, None otherwise
    """
    with p.lock_access_data_base:
        return p.data_base.pop(smiles, None)


# if __name__ == '__main__':
#     load_data_base()
#
#     with open("243_conf/data.json", 'r') as f:
#         data = json.load(f)
#
#     with open("../mcts/configurations/243_conf.json", 'r') as f:
#         config = json.load(f)
#
#     for k, v in data.items():
#         full_smiles = "".join(config["long_prefix"]) + k
#         p.data_base[full_smiles] = v
#
#     save_data_base()
