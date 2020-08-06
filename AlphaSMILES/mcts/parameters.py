from rdkit import RDLogger

# to change RDKit logs
# RDLogger.logger().setLevel(RDLogger.DEBUG)
RDLogger.logger().setLevel(RDLogger.CRITICAL)

f_stop = ""
f_tree = ""
f_info = ""
f_data = ""
f_stat = ""
directory = ""
r_dft = ""

stop = False
models = []
tree = None
tree_info = None
data = None
scorer = None
lock_update_data = None
lock_update_node = None
lock_sa_score = None
data_base = None
lock_access_data_base = None

info_created = "nb_smiles_created"
info_good = "nb_good_smiles"
info_bad = "nb_bad_smiles"
info_alrd_tested = "nb_smiles_already_tested"

id_smile = 0
nb_created = 0
nb_good = 0
turn = 0

# SCORES
s_logp = "logp"
s_sa = "sa"
s_cycle = "cycle"
s_dft = "dft"
s_id = "id"
s_valid = "valid"

tokens = []
config = dict()
