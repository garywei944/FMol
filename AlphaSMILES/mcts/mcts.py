import json
import math
import pickle
import random
import signal
import time
import os
import threading
from joblib import Parallel, delayed

from mcts import parameters as p
from mcts import data_base
from mcts.smiles import SMILES
from mcts.node import Node
from rnn.rnn import load_model


def load_parameters_mcts(config):
    """
    Load previous tree, data and info if the already exist
    otherwise create new ones
    Initialise Scorer and locks for multithreading

    :return: None
    """
    with open('mcts/configurations/' + config + ".json") as c:
        config = json.load(c)
    p.config = config
    p.directory = 'data_out/' + config['configuration_name'] + '/'
    p.f_tree = p.directory + "tree.pckl"
    p.f_info = p.directory + "info.json"
    p.f_data = p.directory + "data.json"
    p.f_stat = p.directory + "stat.csv"
    p.r_dft = p.directory + "dft/"
    p.f_stop = p.directory + "stop.txt"
    p.lock_update_data = threading.Lock()
    p.lock_update_node = threading.Lock()
    p.lock_sa_score = threading.Lock()
    if not os.path.isdir(p.directory) or not os.path.isfile(p.f_tree) or not os.path.isfile(p.f_info)\
            or not os.path.isfile(p.f_stat) or not os.path.isfile(p.f_data):
        if not os.path.isdir(p.directory):
            os.mkdir(p.directory)
        if not os.path.isdir(p.r_dft):
            os.mkdir(p.r_dft)
        with open(p.f_stop, 'w') as stop:
            stop.write("")
        with open(p.f_stat, 'w') as stat:
            stat.write("")
        print("New tree")
        p.tree_info = dict()
        p.tree_info[p.info_created] = 0
        p.tree_info[p.info_good] = 0
        p.tree_info[p.info_bad] = 0
        p.tree_info[p.info_alrd_tested] = 0
        p.data = dict()
        if config['from'] == 'root':
            p.tree = Node(SMILES(config['prefix']))
        else:
            p.tree = Node()
        p.turn = 0

    else:
        print("Loading previous tree, data and info")
        with open(p.f_tree, 'rb') as f:
            pickler = pickle.Unpickler(f)
            p.tree = pickler.load()
            while p.tree.parent:
                p.tree = p.tree.parent
        with open(p.f_data, 'r') as f:
            p.data = json.load(f)
        with open(p.f_info, 'r') as f:
            p.tree_info = json.load(f)
        p.turn = sum(1 for line in open(p.f_stat))

        if config['from'] == 'root':
            if p.tree.smiles.element != config['prefix']:
                raise Exception("Root node different from previous execution")

    with open('rnn_models/' + p.config['rnn_repertory'] + "/config.json") as info_rnn:
        p.tokens = json.load(info_rnn)['tokens']

    data_base.load_data_base()

    reset_score_visit(p.tree)
    p.scorer = getattr(__import__(p.config['scorer'][0], fromlist=[p.config['scorer'][1]]), p.config['scorer'][1])(
        p.config['alpha_scorer'])

    load_scores()

    p.models = load_model(p.config['rnn_repertory'])


def stop_mcts(signum, frame):
    """
    Function called when the process get a SIGINT or SIGTERM
    The MCTS will stop at the end of the current turn

    :param signum: num of the signal
    :param frame: .
    :return: None
    """
    print("The process have been terminated (singal : " + str(signum) + ")")
    print("Please wait for the end of this turn so the current data will be saved")
    print(frame)
    p.stop = True


def signal_handler():
    """
    Handle signal for SIGINT and SIGTERM when the MCTS is running

    :return: None
    """
    signal.signal(signal.SIGINT, stop_mcts)
    signal.signal(signal.SIGTERM, stop_mcts)


def stop_next_turn():
    """
    Dirty way to stop the MCTS in a clean way (without SIGINT or SIGTERM)...
    the mcts finish current turn save data and stop (if you are using dft it can take some time...)
    write "stop" in the file MCTS/stop_mcts

    :return: None
    """
    with open(p.f_stop) as f:
        stop = f.read()
    if "stop" in stop:
        print("MCTS stopped with signal 'stop' in '%s' file" % p.f_stop)
        return True
    return False


def launch():
    """
    Launch the MCTS for nb_turn turns
    It can be stopped in a clean way (saving all current data, tree, ... )
    by writting stop in the file MCTS/stop_mcts
    or sending a SIGINT or SIGTERM
    be careful a SIGKILL kill the process without saving the tree or current data
    if you are using DFT it can takes several hour before stopping it...

    :return: None
    """
    print("Let's begin")
    start = time.time()
    node = p.tree
    node.echo()
    signal_handler()
    prefix = SMILES(p.config['prefix'])
    nb_turn = p.turn + p.config["nb_turn"]
    while (p.turn < nb_turn) and (not stop_next_turn()) and not p.stop:
        last_good = p.tree_info[p.info_good]
        last_bad = p.tree_info[p.info_bad]
        last_already = p.tree_info[p.info_alrd_tested]
        last_created = p.tree_info[p.info_created]
        last_time = time.time()

        print("Turn %d" % p.turn)
        if p.config['from'] == 'node':
            node = get_node_with_prefix(p.tree, prefix)
        node_to_expand = selection(node)
        print("Node to expand : " + str(node_to_expand))
        new_node = expansion(node_to_expand)
        print("New nodes : " + str(new_node))
        new_smiles = simulation(new_node)
        print("new smiles : " + str(new_smiles))
        update(new_smiles)  # (new_node, new_smiles)
        save_tree(node)
        save_data_and_info()
        data_base.save_data_base()

        print("Found %d valid SMILES, %d bad SMILES, %d already tested out of %d generated" % (
            p.tree_info[p.info_good], p.tree_info[p.info_bad], p.tree_info[p.info_alrd_tested],
            p.tree_info[p.info_created]))
        current_created = p.tree_info[p.info_created] - last_created
        current_already = p.tree_info[p.info_alrd_tested] - last_already
        current_good = p.tree_info[p.info_good] - last_good
        current_bad = p.tree_info[p.info_bad] - last_bad
        current_time = int(time.time() - last_time)
        with open(p.f_stat, 'a') as stat:
            stat.write("%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n" % (
                p.turn, p.tree_info[p.info_created], p.tree_info[p.info_good], p.tree_info[p.info_bad],
                p.tree_info[p.info_alrd_tested], current_created, current_good, current_bad, current_already,
                node.get_size(), node.get_height(), int(time.time() - start), current_time))
        p.turn += 1
    print("MCTS worked during %d s" % int(time.time() - start))


def selection_to_delete(node):
    """
    Select the node to use according to the ubc calcul

    :param node: the root node to start the search
    :type node: Node
    :return: the node selected
    """
    while not node.smiles.terminal():
        if not node.children:
            return node
        else:
            node = ubc(node)
    return node


def ubc(node):
    """
    Choose the child depending on the number of visit and the score of the child

    :param node: the parent node
    :type node: Node
    :return: the chosen child
    """
    if not node.children:
        raise Exception("ubc impossible, no child : " + repr(node))
    best_score = -float("inf")
    best_children = []
    for c in node.children:
        if not c.smiles.terminal():
            exploit = c.score / c.visits
            explore = math.sqrt(2.0 * math.log(node.visits) / float(c.visits))
            score = exploit + p.config['exploration_vs_exploitation'] * explore
            if score == best_score:
                best_children.append(c)
            if score > best_score:
                best_children = [c]
                best_score = score
    if not best_children:
        node.visits += 1
        return node.parent
    return random.choice(best_children)


def selection(node):
    """
    Select the node to use according to the ubc calculation

    :param node: the root node to start the search
    :return: the node selected
    """
    while node.children:
        node = ubc(node)
    return node


def expansion(nodes_to_expand):
    """
    Function called for the expansion of the tree,
    Add as much children as the function next_atoms of the current SMILES give

    :param nodes_to_expand: the leaf selected by the function select
    :type nodes_to_expand: Node
    :return: the children created for the leaf
    """
    new_nodes = []
    for s in nodes_to_expand.smiles.next_atoms():
        nodes_to_expand.new_child(s)
        new_nodes.append(nodes_to_expand.children[-1])
    return new_nodes


def simulation(nodes):
    """
    Roll-out,
    Create as much as SMILES for each new created nodes as given in the configuration file

    :param nodes: nodes to test during the simulation part
    :type nodes: list of Node
    :return: smiles created during the simulation
    """
    all_new_smiles = []
    if len(p.models) != 1:
        for n in nodes:
            new_smiles = []
            for i in range(p.config['SMILES_simulated_per_node']):
                smiles = n.smiles.end_smiles_with_model(i)
                while not smiles.terminal():
                    smiles = n.smiles.end_smiles_with_model(i)
                new_smiles.append(smiles)
            all_new_smiles += new_smiles
    else:
        for n in nodes:
            new_smiles = []
            while len(new_smiles) < p.config['SMILES_simulated_per_node']:
                smiles = n.smiles.end_smiles()
                if smiles.terminal():
                    new_smiles.append(smiles)
            all_new_smiles += new_smiles
    return all_new_smiles


def update_smiles(smiles):
    """
    Find the reward for the smiles,
    if its have already been calculated
    then he load it from the data

    :param smiles: the smile to reward
    :type smiles: SMILES
    :return: None
    """
    already = False
    full_smiles = "".join(p.config['long_prefix']) + "".join(smiles.element)
    properties = data_base.select(full_smiles)

    if repr(smiles) in p.data.keys():
        already = True
        with p.lock_update_data:
            p.tree_info[p.info_alrd_tested] += 1
        smiles.properties = p.data[repr(smiles)]
    else:
        if not properties:
            smiles.calculation_of_properties()
            data_base.create(full_smiles, smiles.properties)
        else:
            smiles.properties = properties
        with p.lock_update_data:
            p.data[repr(smiles)] = smiles.properties
    reward = p.scorer.reward(smiles=p.data[repr(smiles)], already=already)
    with p.lock_update_node:
        node = get_node_starting_with("".join(smiles.element))
        # print("Reward %s : %f on node %s" % (str(smiles), reward, repr(node)))
        node.update(reward)


def update(new_smiles):
    """
    Call the function to update nodes with all generated smiles
    Use multitheading if n_jobs in parameters is greater than 1

    :param new_smiles: generated smiles
    :type new_smiles: list of SMILES
    :return: None
    """
    p.tree_info[p.info_created] += len(new_smiles)
    print("Update")
    print("%d smiles to process" % len(new_smiles))
    print("%d unique smiles" % len(set(new_smiles)))
    p.tree_info[p.info_alrd_tested] += (len(new_smiles) - len(set(new_smiles)))
    backend = 'threading'
    Parallel(n_jobs=p.config['n_jobs'], backend=backend)(delayed(update_smiles)(s) for s in set(new_smiles))


def get_node_with_prefix(node, smiles):
    """
    Find or create node to match the smiles
    be carefull, create nodes until the smiles is complete in the tree

    :exemple:

        node = Node()
        node.new_child(SMILES(['O']))
        node.children[0].new_child(SMILES(node.children[0].smiles.element + ['N']))
        node.children[0].new_child(SMILES(node.children[0].smiles.element + ['C']))
        node.children[0].children[1].new_child(SMILES(node.children[0].children[1].smiles.element + ['c']))
        node.children[0].children[1].new_child(SMILES(node.children[0].children[1].smiles.element + ['F']))
        node.children[0].children[1].new_child(SMILES(node.children[0].children[1].smiles.element + ['D']))
        node.children[0].new_child(SMILES(node.children[0].smiles.element + ['c']))
        node.children[0].new_child(SMILES(node.children[0].smiles.element + ['F']))
        node.new_child(SMILES(node.smiles.element + ['C']))
        node.new_child(SMILES(node.smiles.element + ['N']))
        node.new_child(SMILES(node.smiles.element + ['[NH]']))
        pptree.print_tree(node.out_pptree())
        prefix = get_node_with_prefix(node, SMILES(['O', 'C', 'F', 'G']))
        pptree.print_tree(node.out_pptree())
        pptree.print_tree(prefix.out_pptree())

    :param node: the root node to start the prefix
    :type node: Node
    :param smiles: prefix to generate
    :type smiles: SMILES
    :return: finded or created node
    """
    current_node = node
    for i in range(len(smiles.element)):
        next_node = None
        for n in current_node.children:
            if n.smiles.element == smiles.element[:i + 1]:
                next_node = n
                break
        if next_node:
            current_node = next_node
        else:
            current_node.new_child(SMILES(smiles.element[:i + 1]))
            current_node = current_node.children[-1]
    return current_node


def save_tree(node):
    """
    Save the state of the tree

    :param node: the node to save (save all the tree anyway...)
    :type node: Node
    :return: None
    """
    with open(p.f_tree, "wb") as pf:
        pickler = pickle.Pickler(pf)
        pickler.dump(node)
        print("Tree saved")


def load_tree():
    """
    Load the tree

    :return: the root node of the tree
    """
    with open(p.f_tree, 'rb') as pf:
        pickler = pickle.Unpickler(pf)
        node = pickler.load()
    while node.parent:
        node = node.parent
    return node


def save_data_and_info():
    """
    Save the data and the statistics

    :return: None
    """
    with open(p.f_data, 'w') as f:
        json.dump(p.data, f)
    with open(p.f_info, 'w') as f:
        json.dump(p.tree_info, f)


def get_node_starting_with(smiles):
    """
    Find the node in the tree where the SMILES come from

    :param smiles: SMILES to look for
    :type smiles: str
    :return: parent node of the SMILES
    """
    current_node = p.tree
    while current_node.children and not (repr(current_node.smiles)[1:-1] == smiles):
        children = current_node.children
        children.sort(key=lambda x: len(str(x.smiles)), reverse=True)
        modif = False
        for c in children:
            c_smiles = str(c.smiles)[1:-1]
            len_c_smiles = len(c_smiles)
            if smiles[:len_c_smiles] == c_smiles:
                modif = True
                current_node = c
                break
        if not modif:
            return current_node
    return current_node


def load_scores():
    """
    Load the score in the tree for all data already known

    :return: None
    """
    print("Loading scores")
    length_data = len(p.data.keys())
    for i, s in enumerate(p.data.keys()):
        print("\r%d/%d %s" % (i+1, length_data, s), end=' '*(90-len(s)))
        if len(s) > 0 and s[0] != "_":
            node = get_node_starting_with(s)
            node.update(p.scorer.reward(p.data[s], False))
        elif len(s) == 0:
            node = get_node_starting_with(s)
            node.update(p.scorer.reward(p.data[s], False))
    print()


def reset_score_visit(node):
    """
    Reset the scores in the tree

    :param node: root node of the tree
    :type node: Node
    :return: None
    """
    node.visits = 0
    node.score = 0.0
    if node.children:
        for c in node.children:
            reset_score_visit(c)
