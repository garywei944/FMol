import os
import datetime
import json
import numpy as np
import glob
from keras.layers import Dense, TimeDistributed, GRU, Dropout
from keras.layers.embeddings import Embedding
from keras.models import Sequential, model_from_json
from keras.optimizers import Adam
from keras.preprocessing import sequence
from keras.utils.np_utils import to_categorical
from sklearn.utils import resample
from sklearn.base import BaseEstimator, RegressorMixin
from tensorflow.python.keras.callbacks import TensorBoard, EarlyStopping, ModelCheckpoint


def make_filter(tokens_allowed, cn=False):
    """
    This function return a function that return True if all the tokens of the SMILES gived in parameters
    are allowed

    :param tokens_allowed: tokens allowed in your SMILES
    :type tokens_allowed: list of str
    :param cn: if you want C#N and N#C
    :type cn: bool
    :return: function taking a SMILES as parameter and return True if the SMILES is accepted
    """
    if cn:
        def is_smiles_ok(smiles):
            s = smiles
            jump = 0
            for i, c in enumerate(s):
                if jump > 0:
                    jump -= 1
                else:
                    if c not in tokens_allowed:
                        if c == 'N':
                            if i < len(s) - 2:
                                if s[i + 1] == '#' and s[i + 2] == 'C':
                                    jump = 2
                                else:
                                    return False
                            else:
                                return False
                        elif c == 'C':
                            if i < len(s) - 2:
                                if s[i + 1] == '#' and s[i + 2] == 'N':
                                    jump = 2
                                else:
                                    return False
                            else:
                                return False
                        else:
                            return False
            return True
        return is_smiles_ok
    else:
        return lambda smiles: all([x in tokens_allowed for x in smiles])


def parse_data(smiles):
    """
    Convert a list of SMILES into a list of list of tokens and add a '&' at the beginning

    :todo: Replace this function with a parser

    :exemple:

        >>> parse_data(["CCCCn1cccc1", "CCc1c(C)ncnc1O", "CCOC(=C)OCC"])
        [['&', 'C', 'C', 'C', 'C', 'n', '1', 'c', 'c', 'c', 'c', '1', '\\n'],
        ['&', 'C', 'C', 'c', '1', 'c', '(', 'C', ')', 'n', 'c', 'n', 'c', '1', 'O', '\\n'],
        ['&', 'C', 'C', 'O', 'C', '(', '=', 'C', ')', 'O', 'C', 'C', '\\n']]

    :param smiles: the list to convert
    :type smiles: list of str
    :return: list of the converted SMILES
    """
    # add & at the beginning and \n at the end of each SMILES and cut the SMILES into tokens

    all_smiles = []
    element_table = ["C", "N", "B", "O", "P", "S", "F", "Cl", "Br", "I", "(", ")", "=", "#"]

    for i in range(len(smiles)):
        print("\r%d/%d" % (i+1, len(smiles)), end="")
        word_space = smiles[i]
        word = []
        j = 0
        while j < len(word_space):
            word_space1 = []
            if word_space[j] == "[":
                word_space1.append(word_space[j])
                j = j + 1
                while word_space[j] != "]":
                    word_space1.append(word_space[j])
                    j = j + 1
                word_space1.append(word_space[j])
                word_space2 = ''.join(word_space1)
                word.append(word_space2)
                j = j + 1
            else:
                word_space1.append(word_space[j])

                if j + 1 < len(word_space):
                    word_space1.append(word_space[j + 1])
                    word_space2 = ''.join(word_space1)
                else:
                    word_space1.insert(0, word_space[j - 1])
                    word_space2 = ''.join(word_space1)

                if word_space2 not in element_table:
                    word.append(word_space[j])
                    j = j + 1
                else:
                    word.append(word_space2)
                    j = j + 2

        if word[-1] != "\n":
            word.append("\n")
        if word[0] != "&":
            word.insert(0, "&")
        all_smiles.append(list(word))
    print()
    return all_smiles


def find_tokens(all_smiles):
    """
    Find all different tokens from a set of data to get the grammar

    :param all_smiles: list of list of tokens
    :type all_smiles: list of list of tokens
    :return: different tokens in all_smiles
    """
    val = ["\n"]
    for smile in all_smiles:
        for a in smile:
            if a not in val:
                val.append(a)
    return val


def convert_data_to_numbers(tokens, smiles):
    """
    Prepare data for the RNN, convert smiles into arrays of number corresponding to their index in the tokens array

    :param tokens: list of the tokens used
    :type tokens: list of str
    :param smiles: SMILES to convert
    :type smiles: list of str
    :return: input and output of the RNN, input start with a & and and with \\n and output end with two \\n
    """
    x_train = [[tokens.index(atom) for atom in smile] for smile in smiles]
    y_train = [x[1:] + [0] for x in x_train]
    return x_train, y_train


def load_model(config_name):
    """
    Load the architecture and weights of the RNN model(s) in the directory 'config_name'

    :param config_name: name of the configuration (the directory name under rnn_models)
    :type config_name: str
    :return: list of loaded models
    """
    save_directory = 'rnn_models/' + config_name + '/'
    nb_models = len(glob.glob(save_directory + 'model_weights_*'))
    models = []
    for i in range(nb_models):
        model_weights = save_directory + 'model_weights_' + str(i) + '.h5'
        model_architecture = save_directory + 'model_architecture_' + str(i) + '.json'
        if os.path.isfile(model_architecture) and os.path.isfile(model_weights):
            with open(model_architecture, 'r') as f:
                model = model_from_json(f.read())
                # Load weights into the new model
                model.load_weights(model_weights)
                print("Loaded model from disk")
                models.append(model)
        else:
            raise Exception("Model not found at " + save_directory + " " + model_architecture + " " + model_weights)
    return models


def prepare_data(config):
    """
    Load the SMILES file selected in the config file, filter the SMILES with the tokens allowed.
    Then return all selected SMILES and the tokens present.

    :param config: configuration for the RNN
    :type config: dict
    :return: list of the tokens and all selected SMILES
    """
    with open('data_in/' + config["data_input"]) as f_smiles:
        all_smiles = f_smiles.readlines()
    print("Converting SMILES")
    all_smiles = parse_data(all_smiles)
    print("tokens_available = %s" % str(find_tokens(all_smiles)))
    print("Selecting SMILES")
    filter_smiles = make_filter(tokens_allowed=config["tokens_allowed"], cn=config['C#N allowed'])
    all_smiles = list(filter(filter_smiles, all_smiles))
    print("%d SMILES selected" % len(all_smiles))
    tokens = find_tokens(all_smiles)
    config["tokens"] = tokens

    print("tokens length = %d" % len(tokens))
    print("tokens = %s " % str(tokens))
    return tokens, all_smiles


def train_rnn(config, data, tokens, number):
    """
    Train one RNN, keep the best weights of the model and stop it when it doesnt learn anymore

    :param config: config to use
    :type config: dict
    :param data: data to use to train the RNN
    :type data: list of list of str
    :param tokens: list of tokens used to train the RNN
    :type tokens: list of str
    :param number: id of the model
    :type number: int
    :return: None
    """
    print("Model : " + str(number))
    x_train, y_train = convert_data_to_numbers(tokens, data)

    print("SMILES converted to numbers")

    maxlen = 81

    x = sequence.pad_sequences(x_train, maxlen=maxlen, dtype='int32',
                               padding='post', truncating='pre', value=0.)
    y = sequence.pad_sequences(y_train, maxlen=maxlen, dtype='int32',
                               padding='post', truncating='pre', value=0.)

    print("Loading y_train_one_hot")

    y_train_one_hot = np.array([to_categorical(y_i, num_classes=len(tokens)) for y_i in y])
    print(y_train_one_hot.shape)

    n = x.shape[1]

    model = Sequential()

    model.add(Embedding(input_dim=len(tokens), output_dim=len(tokens), input_length=n, mask_zero=False))
    model.add(GRU(units=256, return_sequences=True, activation="tanh", input_shape=(81, 26)))
    model.add(Dropout(0.2))
    model.add(GRU(256, activation='tanh', return_sequences=True))
    model.add(Dropout(0.2))
    model.add(TimeDistributed(Dense(len(tokens), activation='softmax')))
    optimizer = Adam(lr=config['learning_rate'])
    model.summary()

    save_directory = 'rnn_models/' + config['configuration_name'] + '/'
    log = save_directory + 'log'
    model_weights = save_directory + 'model_weights_' + str(number) + '.h5'
    model_architecture = save_directory + 'model_architecture_' + str(number) + '.json'

    if not os.path.isdir(save_directory):
        os.mkdir(save_directory)
    if not os.path.isdir(log):
        os.mkdir(log)

    tensorboard = TensorBoard(log_dir=log + "/{}".format(str(datetime.datetime.now()) + '_model' + str(number)))
    early_stopping = EarlyStopping(monitor='val_loss', patience=10, verbose=0, mode='min')
    mcp_save = ModelCheckpoint(model_weights, save_best_only=True, monitor='val_loss', mode='min')

    model.compile(loss='categorical_crossentropy', optimizer=optimizer, metrics=['accuracy'])

    model.fit(x, y_train_one_hot, epochs=config['epochs'], batch_size=512,
              validation_split=0.1, callbacks=[tensorboard, early_stopping, mcp_save])

    model_json = model.to_json()
    with open(model_architecture, "w") as json_file:
        json_file.write(model_json)
    print("Model saved")

    with open(save_directory + "config.json", 'w') as conf:
        json.dump(config, conf)


def create_rnn(config):
    """
    Load config and create the RNN

    :param config: name of the config to use
    :type config: str
    :return: None
    """
    with open('rnn/configurations/' + config + '.json') as c:
        config = json.load(c)

    tokens, all_smiles = prepare_data(config)

    if config['bootstrapping'] < 1:
        train_rnn(config, all_smiles, tokens, 0)
    else:
        boots = []
        for i in range(config['bootstrapping']):
            boots.append(resample(all_smiles, replace=False, n_samples=config['nb_samples_rnn_bt'], random_state=i))

        for i in range(config['bootstrapping']):
            train_rnn(config, boots[i], tokens, i)


class Model(BaseEstimator, RegressorMixin):
    """
    Model using Keras
    """
    def __init__(self, lr=0.1, epochs=300, tokens=None, n=0):
        """
        Init of the model layers

        :param lr: learning rate
        :type lr: float
        :param epochs: number of epochs
        :type epochs: int
        :param tokens: tokens used to train the RNN
        :type tokens: list of str
        :param n: input_length
        :type n: int
        """
        if tokens is None:
            tokens = ['\n', '&', 'C', '1', '2', '3', '4', 'O', 'N', '(', ')', '=', 'c', '[nH]', 'S', 'n', 's', 'o', '#', 'Cl', '[NH]']
        self.lr = lr
        self.epochs = epochs
        self.number = datetime.datetime.now()
        self.estimator = Sequential()
        self.estimator.add(Embedding(input_dim=len(tokens), output_dim=len(tokens), input_length=n, mask_zero=False))
        self.estimator.add(GRU(units=256, return_sequences=True, activation="tanh", input_shape=(81, 26)))
        self.estimator.add(Dropout(0.2))
        self.estimator.add(GRU(256, activation='tanh', return_sequences=True))
        self.estimator.add(Dropout(0.2))
        self.estimator.add(TimeDistributed(Dense(len(tokens), activation='softmax')))
        self.tensorboard = None
        self.early_stopping = None
        self.mcp_save = None

    def fit(self, x, y):
        """
        Compile and fit the model

        :param x: training data
        :type x: list of list of int
        :param y: target data
        :type y: list of list of int
        :return: None
        """
        optimizer = Adam(lr=self.lr)
        self.estimator.summary()

        save_directory = 'rnn_models/' + "test" + '/'
        log = save_directory + 'log'
        model_weights = save_directory + 'model_weights_' + str(self.number) + '.h5'
        # model_architecture = save_directory + 'model_architecture_' + str(number) + '.json'

        if not os.path.isdir(save_directory):
            os.mkdir(save_directory)
        if not os.path.isdir(log):
            os.mkdir(log)

        self.tensorboard = TensorBoard(
            log_dir=log + "/{}".format(str(datetime.datetime.now()) + '_model' + str(self.number)))
        self.early_stopping = EarlyStopping(monitor='val_loss', patience=10, verbose=0, mode='min')
        self.mcp_save = ModelCheckpoint(model_weights, save_best_only=True, monitor='val_loss', mode='min')

        self.estimator.compile(loss='categorical_crossentropy', optimizer=optimizer, metrics=['accuracy'])
        self.estimator.fit(x, y, epochs=self.epochs, batch_size=512, validation_split=0.1,
                           callbacks=[self.tensorboard, self.early_stopping, self.mcp_save])

    def score(self, x, y, sample_weight=None):
        """
        Score the model

        :param x: training data
        :type x: list of list of int
        :param y: target data
        :type y: list of list of int
        :param sample_weight:
        :return: Score for the model
        """
        score = self.estimator.evaluate(x, y)
        print(score)
        return score[1]
