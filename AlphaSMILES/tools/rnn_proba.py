import json
import os

import numpy as np
import tensorflow as tf
from keras.preprocessing import sequence

from mcts import parameters as p
from mcts.smiles import mol_to_int
from rnn.rnn import load_model

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
# tf.compat.v1.logging.set_verbosity(50)
tf.logging.set_verbosity(tf.logging.ERROR)
os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
os.environ["CUDA_VISIBLE_DEVICES"] = ""

os.chdir('..')


def proba_different_position(model, smiles_to_test, index, atom, verbose=False):
    """
    Print the probability of apparition of a SMILES with different position of an atom

    :exemple:

        >>> proba_different_position(models[0],
        >>>                          ['c', '1', 'c', 'c', 'c', 'c', '1', '\\n'],
        >>>                          [0, 2, 3, 4, 5], 's', False)
        s1cccc1 : 0.00000000045901715348268789499561297600283699127388103
        c1sccc1 : 0.00000000064832566091357019325257395698722251431433960
        c1cscc1 : 0.00000001438888019436213103261810813092436145410601966
        c1ccsc1 : 0.00000065275450335318335220226868345627657674867805326
        c1cccs1 : 0.00000011832597629642405283315926981316246191511254437

    :param model: model to use
    :type model: Keras Model
    :param smiles_to_test: list of SMILES to test
    :type smiles_to_test: list of str
    :param index: index where the atom heve to be placed
    :type index: list of int
    :param atom: atom to place
    :type atom: str
    :param verbose: True : print all information of each turn, False : print only the probability
    :type verbose: bool
    :return: None
    """
    smiles = []
    for l in index:
        e = smiles_to_test.copy()
        e[l] = atom
        smiles.append(e)
    proba_different_smiles(model, smiles, verbose)


def proba_different_smiles(model, smiles_to_test, verbose=False):
    """
    Print the probability of apparition of the SMILES with the model given in parameter

    :param model: model to use
    :type model: Keras Model
    :param smiles_to_test: list of SMILES to test
    :type smiles_to_test: List[List[[str]]]
    :param verbose: True : print all information of each turn, False : print only the probability
    :type verbose: bool
    :return: None
    """
    for e in smiles_to_test:
        proba = 1
        if verbose:
            print("".join(e))
        for i, l in enumerate(e):
            smiles = ['&'] + e[:i]
            smiles_int = mol_to_int(smiles)
            x = np.reshape(smiles_int, (1, len(smiles_int)))
            x_pad = sequence.pad_sequences(x, maxlen=81, dtype='int32', padding='post', truncating='pre', value=0.)
            predictions = model.predict(x_pad)
            preds = np.asarray(predictions[0][len(smiles_int) - 1]).astype('float64')
            index_suiv = p.tokens.index(e[i])
            proba = proba * preds[index_suiv]
            if verbose:
                print("SMILES : %s" % "".join(e[:i]))
                print("Suivant : %s" % p.tokens[index_suiv])
                print("Proba suivant : %.9f" % proba)
                for j, pred in enumerate(preds):
                    print("%s\t%.9f" % (p.tokens[j], pred))
                print("Sum preds : %.12f" % np.sum(preds))
        if e[-1] == '\n':
            print("%30s : %.100f" % ("".join(e[:-1]), proba))
        else:
            print("%30s : %.100f" % ("".join(e), proba))


if __name__ == '__main__':
    config_name = 'rnn_2.1'

    smiles = [['c', '1', 'c', 'c', 'c', '(', 's', '1', ')', 'c', '1', 'c', 'c', 'c', 's', '1', '\n'],
              ['c', '1', 'c', 'c', '(', 's', 'c', '1', ')', 'c', '1', 's', 'c', 'c', '(', 'c', '1', ')', '\n'],
              ['c', '1', 'c', 'c', 'c', '(', 's', '1', ')', 'c', '1', 's', 'c', 'c', '(', 'c', '1', ')', '\n'],
              ['c', '1', 'c', 'c', '(', 's', 'c', '1', ')', 'c', '1', 'c', 'c', 'c', 's', '1', '\n'],
              ['s', '1', 'c', 'c', 'c', 'c', '1', 'c', '2', 'c', 'c', 'c', 's', '2', '\n']
              ]

    models = load_model(config_name)

    with open('rnn_models/' + config_name + "/config.json") as info_rnn:
        p.tokens = json.load(info_rnn)['tokens']

    print(p.tokens)

    proba_different_position(models[0], ['c', '1', 'c', 'c', 'c', 'c', '1', '\n'], [0, 2, 3, 4, 5], 's', False)

    # proba_different_smiles(models[0], smiles, False)

    # from rnn.rnn import convert_data_to_numbers
    # from keras.utils.np_utils import to_categorical
    #
    # x_train, y_train = convert_data_to_numbers(p.tokens, [['&', 'C', '\n']])
    #
    # maxlen = 81
    #
    # x = sequence.pad_sequences(x_train, maxlen=maxlen, dtype='int32',
    #                            padding='post', truncating='pre', value=0.)
    # y = sequence.pad_sequences(y_train, maxlen=maxlen, dtype='int32',
    #                            padding='post', truncating='pre', value=0.)
    #
    # print("Loading y_train_one_hot")
    #
    # y_train_one_hot = np.array([to_categorical(y_i, num_classes=len(p.tokens)) for y_i in y])
    # print(y_train_one_hot.shape)
    #
    # from keras.optimizers import Adam
    # models[0].compile(optimizer=Adam(lr=0.001), loss='categorical_crossentropy')
    #
    # for i in range(0,10):
    #     proba_different_smiles(models[0], [['&', 'C', '\n']], False)
    #     print(models[0].evaluate(x, y_train_one_hot))


