Look for the best RNN parameters
--------------------------------

You can execute a grid search to find the best parameters for the RNN. This script allow you to find a good learning rate for your RNN but you can also implement other test by changing the *Model* class in the *'rnn/rnn.py'* file to get access to other parameters.



.. code-block:: python
    :linenos:

    import json
    import numpy as np
    from keras.preprocessing import sequence
    from keras.utils.np_utils import to_categorical
    from rnn.rnn import prepare_data

    with open('rnn/configurations/rnn_test.json') as c:
        config = json.load(c)

    tokens, all_smiles = prepare_data(config)
    #%%
    from rnn.rnn import convert_data_to_numbers
    x_train, y_train = convert_data_to_numbers(tokens, all_smiles[:100000])

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
    #%%
    from rnn.rnn import Model
    from sklearn.model_selection import GridSearchCV
    model = Model(tokens=tokens, n= x.shape[1])
    grid = {'lr': [0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001]}
    gridsearch = GridSearchCV(model, grid, cv=10)
    gridsearch.fit(x, y_train_one_hot)
    #%%
    gridsearch.best_params_


Launching this script in a jupyter-notebook is maybe easier for analysing results.
