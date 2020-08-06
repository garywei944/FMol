import hashlib
import numpy as np
from keras.preprocessing import sequence

from mcts import parameters as p


class SMILES:
    """
    Class representing a SMILES
    """

    def __init__(self, element=None):
        """
        Initialise the new SMILES

        :param element: the list of token corresponding to the new SMILES
        :type element: list of str
        """
        self.element = element if element else []
        if element and element[-1] == '\n':
            self.properties = dict()
            self.scores = dict()
        else:
            self.properties = dict()
            self.scores = None

    def next_atoms(self):
        """
        Generate new SMILES with one more token than the current SMILES by using the RNN model
        SMILES are selected is their probability is greater than 0.0001% according to the RNN model

        :return: the list of more probables next SMILES
        """
        smiles_to_expand = []
        if p.config['expansion'] == "all":
            for v in p.tokens[2:]:
                if v == ")" and (self.element.count("(") <= self.element.count(')')):
                    pass
                elif v in ["1", "2", "3", "4", "=", "#", "("] and not self.element:
                    pass
                elif v == "2" and (self.element.count("1") == 0):
                    pass
                elif v == "3" and ((self.element.count("1") == 0) and (self.element.count("2") == 0)):
                    pass
                elif v == "4" and ((self.element.count("1") == 0) and (self.element.count("2") == 0)
                                   and (self.element.count("3") == 0)):
                    pass
                elif v in ["#", "="] and (self.element[-1] in ["#", "=", 'Cl']):
                    pass
                elif v == "Cl" and self.element and (self.element[-1] in ["#", "="]):
                    pass
                else:
                    smiles_to_expand.append(SMILES(self.element + [v]))
        elif p.config['expansion'] == "best":
            smiles_int = mol_to_int(['&'] + self.element)
            x = np.reshape(smiles_int, (1, len(smiles_int)))
            x_pad = sequence.pad_sequences(x, maxlen=81, dtype='int32',
                                           padding='post', truncating='pre', value=0.)
            predictions = p.models[0].predict(x_pad)
            preds = np.asarray(predictions[0][len(smiles_int) - 1]).astype('float64')
            preds = preds / np.sum(preds)
            next_probas = np.random.multinomial(1, preds[1:], 30)
            next_int = np.argmax(next_probas, axis=1) + 1
            to_expand = list(set(next_int))
            for s in to_expand:
                smiles_to_expand.append(SMILES(int_to_smile(smiles_int[1:]) + [p.tokens[s]]))
        elif p.config['expansion'] == "proba":
            smiles_int = mol_to_int(['&'] + self.element)
            x = np.reshape(smiles_int, (1, len(smiles_int)))
            x_pad = sequence.pad_sequences(x, maxlen=81, dtype='int32',
                                           padding='post', truncating='pre', value=0.)
            predictions = p.models[0].predict(x_pad)
            preds = np.asarray(predictions[0][len(smiles_int) - 1]).astype('float64')
            smiles_to_expand = []
            for i, pred in enumerate(preds):
                '''
                for the first node:     with 0.0001 10/21 not accepted
                \n      0.000069        no
                &       0.000000        no
                C       0.844086
                O       0.080062
                (       0.000476
                =       0.000779
                )       0.000053        no
                c       0.010339
                1       0.000086        no
                N       0.058904
                n       0.000297
                2       0.000054        no
                3       0.000015        no
                4       0.000002        no
                [nH]    0.000020        no
                Cl      0.001153
                S       0.000430
                o       0.000019        no
                #       0.000123
                [NH]    0.003027
                s       0.000006        no
                '''
                if pred > p.config['proba_min'] and p.tokens[i] != "&" and p.tokens[i] != "\n":
                    smiles_to_expand.append(SMILES(int_to_smile(smiles_int[1:]) + [p.tokens[i]]))
            if not smiles_to_expand:
                for i, pred in enumerate(preds):
                    if pred > (p.config['proba_min'] / 100) and p.tokens[i] != "&" and p.tokens[i] != "\n":
                        smiles_to_expand.append(SMILES(int_to_smile(smiles_int[1:]) + [p.tokens[i]]))
        return smiles_to_expand

    def next_atom(self):
        """
        Generate one new SMILES with one more token than the current SMILES by using the RNN model

        :return: the next SMILES
        """
        smiles_int = mol_to_int(['&'] + self.element)
        x = np.reshape(smiles_int, (1, len(smiles_int)))
        x_pad = sequence.pad_sequences(x, maxlen=81, dtype='int32',
                                       padding='post', truncating='pre', value=0.)
        predictions = p.models[0].predict(x_pad)
        preds = np.asarray(predictions[0][len(smiles_int) - 1]).astype('float64')

        preds = preds[2:]
        preds = preds / np.sum(preds)
        next_probas = np.random.multinomial(1, preds, 1)
        next_int = np.argmax(next_probas) + 2
        smiles_int.append(next_int)
        return SMILES(int_to_smile(smiles_int[1:]))

    def end_smiles(self):
        """
        Complete the current SMILES by using the RNN model
        Add tokens until getting a '\\n'

        :return: the new SMILES generated
        """
        smiles_int = mol_to_int(['&'] + self.element)
        while not smiles_int[-1] == p.tokens.index("\n"):
            x = np.reshape(smiles_int, (1, len(smiles_int)))
            x_pad = sequence.pad_sequences(x, maxlen=81, dtype='int32',
                                           padding='post', truncating='pre', value=0.)
            predictions = p.models[0].predict(x_pad)
            preds = np.asarray(predictions[0][len(smiles_int) - 1]).astype('float64')
            preds = preds / np.sum(preds)
            next_probas = np.random.multinomial(1, preds, 1)
            next_int = np.argmax(next_probas)
            smiles_int.append(next_int)
            if len(smiles_int) > 81:
                break
        return SMILES(int_to_smile(smiles_int[1:]))

    def end_smiles_with_model(self, index_model):
        """
        Complete the SMILES until reaching '\\n' with the model given

        :param index_model: index of the model
        :type index_model: int
        :return: the new SMILES generated
        """
        smiles_int = mol_to_int(['&'] + self.element)
        while not smiles_int[-1] == p.tokens.index("\n"):
            x = np.reshape(smiles_int, (1, len(smiles_int)))
            x_pad = sequence.pad_sequences(x, maxlen=81, dtype='int32',
                                           padding='post', truncating='pre', value=0.)
            predictions = p.models[index_model].predict(x_pad)
            preds = np.asarray(predictions[0][len(smiles_int) - 1]).astype('float64')
            preds = preds / np.sum(preds)
            next_probas = np.random.multinomial(1, preds, 1)
            next_int = np.argmax(next_probas)
            smiles_int.append(next_int)
            if len(smiles_int) > 81:
                break
        return SMILES(int_to_smile(smiles_int[1:]))

    def terminal(self):
        """
        Return True if the current SMILES is complete (if there is a '\\n' at the end of the SMILES)

        :return: True if the SMILES is terminal
        """
        return (self.element[-1] == '\n') if self.element else False

    def calculation_of_properties(self):
        """
        Calculation of different properties asked in the config file.

        :return: None
        """
        from mcts.properties.properties import Property2D
        prop = Property2D(self)
        for imp, pro in p.config['properties']:
            prop = prop << getattr(__import__(imp, fromlist=[pro]), pro)
        prop.calculate()

    def __eq__(self, other):
        """
        Equality between two SMILES

        :param other: the other SMILES
        :type other: SMILES
        :return: boolean if the two SMILES are the same or not
        """
        return repr(self) == repr(other)

    def __hash__(self):
        """
        Hash the SMILES

        :return: the hash of the SMILES
        """
        return hash(hashlib.sha256("".join(self.element).encode('utf-8')).hexdigest())

    def __repr__(self):
        """
        Representation of the SMILES
        The SMILES is surrounded with _ if its not terminal

        :return: a string representing the SMILES
        """
        if self.element and self.element[-1] == '\n':
            return "".join(self.element[:-1])
        else:
            return "_" + "".join(self.element) + "_"


def int_to_smile(list_of_int):
    """
    Convert the list of int to list of token according to the vocabulary

    :param list_of_int: list of int representing the SMILES
    :type list_of_int: list of int
    :return: list of token
    """
    return [p.tokens[s] for s in list_of_int]


def mol_to_int(list_of_mol):
    """
    Convert the list of tokens to list of int according to their index in the vocabulary

    :param list_of_mol: list of tokens
    :type list_of_mol: list of str
    :return: list of indices of the tokens
    """
    return [p.tokens.index(s) for s in list_of_mol]
