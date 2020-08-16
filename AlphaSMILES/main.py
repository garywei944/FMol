import tensorflow as tf

no_gpu = False
train_rnn = True
launch_mcts = True
no_warning = True

if no_warning:
    import os

    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
    # tf.logging.set_verbosity(tf.logging.ERROR)

if no_gpu:
    import os

    os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
    os.environ["CUDA_VISIBLE_DEVICES"] = ""
else:
    physical_devices = tf.config.list_physical_devices('GPU')
    tf.config.experimental.set_memory_growth(physical_devices[0], enable=True)

if train_rnn:
    from rnn.rnn import create_rnn

    create_rnn('rnn_sample')

if launch_mcts:
    from mcts.mcts import load_parameters_mcts, launch

    load_parameters_mcts('mcts_sample')
    launch()
