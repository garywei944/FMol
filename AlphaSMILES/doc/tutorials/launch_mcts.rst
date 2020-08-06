Launch your MCTS
================

Configuration
-------------

configuration_name
::::::::::::::::::

The configuration name will be the name of the output data under the directory *'data_out/'*

.. code-block:: python

    config['configuration_name'] = "mcts_test"

rnn_repertory
:::::::::::::

This is the name of the rnn configuration (the name of the RNN directory under *'rnn_models/'*).

.. code-block:: python

    config['rnn_repertory'] = "rnn_test"

long_prefix
:::::::::::

The long prefix is a prefix the the MCTS and the RNN won't see, only the calculation of properties will add it. You should use it when your prefix is more than 15-20 tokens. The prefix have to be a complete SMILES without the '\\n'.

.. code-block:: python

    config['long_prefix'] = ['c', '1', 'c', '2', 'c', '(', '=', 'O', ')', 'n', '(', 'C', ')', 'c',
                             '(', '=', 'O', ')', 'c', '(', 'c', 'c', '3', ')', 'c', '2', 'c', '4',
                             'c', '3', 'c', '2', 'c', 'c', 'c', 'c', 'c', '2', 's', 'c', '4', 'c',
                             '1']

prefix
::::::

There is two ways of managing short prefix.

+ from root, the prefix will be the prefix of the tree and you will be able to use the generated tree only with this prefix.

+ from node, the tree will be expand to your prefix (or will select the prefix if it already exist in the tree) and will always use this node to start.

Represent the prefix like this : ['c', '1']

Note that all tokens present in the prefix have to be in the RNN tokens.

.. code-block:: python

    config['prefix'] = []
    config['from'] = 'node'  # 'node' or 'root'


SMILES_simulated_per_node
:::::::::::::::::::::::::

This parameter give the number of SMILES that will be created for each turn of the MCTS.

If you are using bootstrapping, if the number is inferior at the number of model, only the *n* first one will be used, if the number is greater the firsts model could be reused.

Keep in mind that using a big value here won't increase the yield of the MCTS, on the contrary, the MCTS will take more time to converge to a set of good solutions.

.. code-block:: python

    config['SMILES_simulated_per_node'] = 2


nb_turn
:::::::

Number of turn that the MCTS have to do before stopping. You can also tell him to stop at the end off his turn by sending it a SIGINT or SIGTERM (if you are using DFT calculation, writing *stop* in the file *'data_out/[config_name]/stop.txt'* will stop the MCTS at the end of the turn too).

.. code-block:: python

    config['nb_turn'] = 200

exploration_vs_exploitation
:::::::::::::::::::::::::::

The number managing the ratio between exploration and exploitation help to encourage the MCTS to explore new branch if the number is higher. If its lower, it will maximise the exploitation of a branch that gives good scores.

.. code-block:: python

    config['exploration_vs_exploitation'] = 1

expansion and proba_min
:::::::::::::::::::::::

Here you can choose the expansio, you have 3 choices:

+ *"proba"* : the MCTS will use the probability given by the RNN model and will open the nodes with the tokens with a probability higher than *"proba_min"*.

+ *"all"* : the MCTS will open all nodes according to rules, those rules are incomplete so a lot of nodes are open on bad solution (impossibles SMILES)

+ *"best"* : the MCTS will open the best choice given by the RNN model (sometimes only 2 or 3 nodes).

The systems *"best"* and *"proba"* rely too much on the RNN model, so when the RNN is in difficulty (generating a SMILES from a big prefix) the RNN make a lot of mistakes.

The system *"all"* is currently unfinished and open to many bad nodes.

.. code-block:: python

    config['expansion'] = "proba"
    config['proba_min'] = 0.0001

n_jobs and nb_core_dft
::::::::::::::::::::::

This two parameters help you to parallelize the execution.

*'n_jobs'* is the number of thread created during the update( when the properties of the SMILES are calculated).

*'nb_core_dft'* is the number of core used during the DFT calculation with Gaussian.

.. code-block:: python

    # multitasking
    config['n_jobs'] = 10
    config['nb_core_dft'] = 8

properties
::::::::::

The properties are calculated during the update to help you give a score to a molecule. See `create your properties <create_your_properties.html>`_ to create new ones.

This have to be a tuple, the first entry is the location of the property, the second is the name of the class.

The property checking if a SMILES is good or bad will always be checked.

.. code-block:: python

    config['properties'] = [("mcts.properties.properties", "SAScoreProperty2DDecorator"),
                            ("mcts.properties.properties", "CycleProperty2DDecorator"),
                            ("mcts.properties.properties", "LogPProperty2DDecorator"),
                            ("mcts.properties.properties", "DFTPropertyDecorator"),
                            ]

scorer and alpha_scorer
:::::::::::::::::::::::

The scorer use the properties of a SMILES to give it a grade.

The *'alpha_scorer'* if the *alpha* of the function how calculate the reward from the score (between 1 and -1).

.. math::

    reward = { | alpha * score | \over 1 + | alpha * score |}

See `create your scorer <create_your_scorer.html>`_ to learn how to create a new scorer.

.. code-block:: python

    config['scorer'] = ("mcts.scorer.scorer", 'ScorerValidSMILES')
    config['alpha_scorer'] = 1

data_base
:::::::::

This is the name of the data base with the results of all generated molecules located in the folder *'data_out/'*

.. code-block:: python

    config["data_base"] = "3-21G"

Launch the MCTS
---------------

Once you have generated a json configuration file, you can use the file *main.py* to launch the MCTS with the desire configuration or use that kind of script:


.. code-block:: python
    :linenos:
    :emphasize-lines: 13,14,15

    # To disable TensorFlow warnings
    import os
    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
    import tensorflow as tf
    tf.logging.set_verbosity(tf.logging.ERROR)

    # If you want to work with the GPU, set no_gpu to False
    no_gpu = True
    if no_gpu:
        os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
        os.environ["CUDA_VISIBLE_DEVICES"] = ""

    from mcts.mcts import load_parameters_mcts, launch
    load_parameters_mcts('mcts_test')
    launch()

After executing this code, generated data will be keep in the directory *'data_out/[configuration name]/'*.
This directory contain :

+r*'dft/'* : a directory for DFT logs

+ *'data.json'* : json file containing all generated SMILES and their properties.

+ *'info.json'* : json file containing the number of SMILES generated, good ones, bad ones and already created ones.

+ *'stat.csv'* : csv file containing statistics about SMILES generated at each turn, time, etc...

+ *'stop.txt'* : text file, if you write *stop* in this file the MCTS will stop at the end of the current turn.

+ *'tree.pckl'* : pickle file containing the tree.
