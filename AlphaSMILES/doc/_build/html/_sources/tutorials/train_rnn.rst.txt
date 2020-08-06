Train your RNN
==============


Configure the RNN
-----------------

In order to train the RNN you need to create a JSON file contaning information about the parameters. You can find a tool to help you build it in the folder *tools/*, *config_creator_rnn.py*.

After setting all parameters execute the script to create a json file with those parameters in *'rnn/configuration'*.

configuration_name
::::::::::::::::::

The configuration name will be the name of the folder containing the model. You will reuse it while configuring the MCTS or when looking for the RNN.

.. code-block:: python

	config['configuration_name'] = "name"


data_input
::::::::::

Data input is the name of the file containing the SMILES to use to train the RNN. This file have to be located in the folder *data_in*. You can find the differents files `here <prepare_your_data.html>`_.

.. code-block:: python

	config['data_input'] = "all_smiles_unique"


tokens_allowed
::::::::::::::

Tokens allowed is the list of tokens that the RNN will use to filter the SMILES file and select the SMILES used to train the RNN.

.. code-block:: python

	config['tokens_allowed'] = ['\n', '&', 'C', 'O', '(', '=', ')', 'c', '1', 'N', 'n', 
				    '2', '3', '4', '[nH]', 'Cl', 'S', 'o', '#', '[NH]', 's']

Here is a list of tokens existing in the *all_smiles_unique* file

.. code-block:: python

	['\n', '&', 'C', 'O', '(', '=', ')', 'c', '1', 'N', 'S', '[O]', '2', '-', 'n', '3', '4', 's',
	'[N+]', '[O-]', 'F', '[Si]', '[nH]', 'Cl', 'Br', '[C]', '[N]', '[PH]', '[S]', 'P', 'o', 'B',
	'[SH]', '#', '[CH]', '5', '6', '[As]', '[NH]', '[Ar]', '[CH2]', '[Se]', '[PH2]', '[AlH3]',
	'[GeH4]', '7', '[SiH]', '[SeH]', '[SeH2]', '[CH3]', '[Cl]', '[Ge]', '[AsH3]', '[Br]',
	'[AsH2]', '[GeH]', '[P]', '[c]', '[Ne]', '[SiH3]', '[Mg]', '[AsH]', '[se]', '[Na]', '[OH]',
	'[SiH4]', '[SiH2]', 'p', '[F]', '[Li]', '[GeH2]', '[Be]', '[GeH3]', '[siH]', '[PH3]', '[si]',
	'[SH3]', '[pH]', '[SH2]', '[NH2]', '[B]', '[PH4]', '8']

When using the *prepare_data* function you can see the list printed at the end of the execution.


.. autofunction:: rnn.rnn.prepare_data
    :noindex:


C#N allowed
:::::::::::

This parameter is used to add 'C#N' and 'N#C' to the possibilities, if 'C', '#' and 'N' are not allowed but you still want 'C#N' in the possibilities.

.. code-block:: python

	config['C#N allowed'] = False

bootstrapping and nb_samples_rnn_bt
:::::::::::::::::::::::::::::::::::

You have the possibility to train and use multiples RNN.
Bootstrapping consist of trainning multiples RNN with differents data by picking
*'nb_samples_rnn_bt'* SMILES in the SMILES file.

The parameter *bootstrapping* allow you to choose if you want to use bootstrapping

.. code-block:: python

	# number of rnn to train
	# 0 : 1 classic rnn trained with all data
	# 1 : 1 rnn trained with 'nb_samples_rnn_bt' samples
	# x : x rnn trained with 'nb_samples_rnn_bt' samples
	config['bootstrapping'] = 0
	config['nb_samples_rnn_bt'] = 300

epochs
::::::

Choose the number max of epochs for the RNN. The RNN will stop before this number when he doesn't learn anymore and will save the best RNN.

.. code-block:: python

	# number of epochs max used to train the rnn
	config['epochs'] = 200

learning_rate
:::::::::::::

Choose the learning rate for the RNN, after some test we conclude that the one that give more regular RNN is 0.001. (see ` how to do a grid search <hyperparameter_rnn.html>`_)

.. code-block:: python

	# learning rate used to train the rnn
	config['learning_rate'] = 0.01


Train it
--------

Once you have generated a json configuration file, you can use the file *main.py* to train the rnn with the desire configuration or use that kind of script:

.. code-block:: python
    :linenos:
    :emphasize-lines: 13,14


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

    from rnn.rnn import create_rnn
    create_rnn('rnn_2.1')

This code disable the warnings from TensorFlow and the use of GPU then create a RNN from the configuration *'rnn_2.1'*.

After executing this code, you will see the progression of the SMILES selection, then will be printed the number of selected SMILES and the tokens presents.
The RNN will be trained until its not learning anything anymore, the best weights are saved as the architecture in the repertory under *'rnn_models/'*
with the same name as the configuration name.

You can watch the performance of the RNN by looking at the graphs generated with TensorsBoard by using the following line (while being in the repertory *'rnn_models/config_name'*
, don't forget to use the correct conda environment):

.. code-block:: bash

    python -m tensorboard.main --logdir=`pwd`/log

