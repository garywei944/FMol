import os

# Set if you want to run AlphaSMILES or AlphaFold
run_alphasmiles = True
run_alphafold = True

# Please set up the following variables according to your PGU
gpu = True
gpu_more_than_8g_memory = False

# AlphaSMILES Configuration
# If you run this project for the first time ignore the following variables. For more specific details, check readme
train_rnn = True
launch_mcts = True
no_warning = True

# AlphaFold Configuration
# Please unzip the weight folders 873731, 916425, and 941521 into 'alphafold_pytorch/model' folder.
# The files can be downloaded at https://storage.googleapis.com/alphafold_casp13_data/alphafold_casp13_weights.zip
target = "T1019s2"
target_file = "test_data/{}.pkl".format(target)
model_dir = "model"

if run_alphasmiles:
    if os.fork():
        os.wait()
    else:
        # from AlphaSMILES.main import launch_alphasmiles
        # os.chdir("AlphaSMILES")
        # launch_alphasmiles(no_gpu=not gpu, train_rnn=train_rnn, launch_mcts=launch_mcts, no_warning=no_warning)
        print("haha")
        exit()

if run_alphafold:
    if os.fork():
        os.wait()
    else:
        from alphafold_pytorch.main import launch_alphafold
        os.chdir("alphafold_pytorch")
        launch_alphafold(target=target, target_file=target_file, model_dir=model_dir,
                         gpu_more_than_8g_memory=gpu_more_than_8g_memory)
        exit()
