no_gpu = False
train_rnn = True
launch_mcts = False
no_warning = True

if no_warning:
    import os

    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
    import tensorflow as tf
    # tf.logging.set_verbosity(tf.logging.ERROR)

if no_gpu:
    import os

    os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
    os.environ["CUDA_VISIBLE_DEVICES"] = ""

if train_rnn:
    from rnn.rnn import create_rnn

    create_rnn('rnn_sample')

if launch_mcts:
    from mcts.mcts import load_parameters_mcts, launch

    load_parameters_mcts('mcts_test')
    launch()
