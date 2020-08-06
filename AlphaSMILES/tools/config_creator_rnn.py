import json

'''
This code help you to prepare the json file of configuration to train the rnn
'''

config = dict()

# name of your configuration
config['configuration_name'] = "rnn_sample"

# name of the file containing all SMILES in the folder 'data_in'
# one SMILES per line
config['data_input'] = "all_smiles"
# config['data_input'] = "1k_rndm_zinc_drugs_clean.smi"

'''
those tokens are the tokens available with the file 'data_in/all_smiles'
tokens_available = ['\n', '&', 'C', 'O', '(', '=', ')', 'c', '1', 'N', 'S', '[O]', '2', '-', 'n', '3', '4', 's',
                    '[N+]', '[O-]', 'F', '[Si]', '[nH]', 'Cl', 'Br', '[C]', '[N]', '[PH]', '[S]', 'P', 'o', 'B',
                    '[SH]', '#', '[CH]', '5', '6', '[As]', '[NH]', '[Ar]', '[CH2]', '[Se]', '[PH2]', '[AlH3]',
                    '[GeH4]', '7', '[SiH]', '[SeH]', '[SeH2]', '[CH3]', '[Cl]', '[Ge]', '[AsH3]', '[Br]',
                    '[AsH2]', '[GeH]', '[P]', '[c]', '[Ne]', '[SiH3]', '[Mg]', '[AsH]', '[se]', '[Na]', '[OH]',
                    '[SiH4]', '[SiH2]', 'p', '[F]', '[Li]', '[GeH2]', '[Be]', '[GeH3]', '[siH]', '[PH3]', '[si]',
                    '[SH3]', '[pH]', '[SH2]', '[NH2]', '[B]', '[PH4]', '8']
'''

# choose the tokens to use to train the rnn
config['tokens_allowed'] = ["\n", "&", "C", "O", "(", "=", ")", "c", "1", "N", "n",
                            "2", "3", "4", "[nH]", "Cl", "S", "o", "#", "[NH]", "s"]

# this is used to add 'C#N' and 'N#C' to the possibilities
config['C#N allowed'] = False

# number of rnn to train
# 0 : 1 classic rnn trained with all data
# 1 : 1 rnn trained with 'nb_samples_rnn_bt' samples
# x : x rnn trained with 'nb_samples_rnn_bt' samples
config['bootstrapping'] = 1
config['nb_samples_rnn_bt'] = 200000

# number of epochs max used to train the rnn
config['epochs'] = 2

# learning rate used to train the rnn
config['learning_rate'] = 0.001

with open('../rnn/configurations/' + config['configuration_name'] + '.json', 'w') as conf:
    json.dump(config, conf)
