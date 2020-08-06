import pptree

from mcts import parameters as p
from mcts.smiles import SMILES


class Node:
    """
    Representation of a node in the MCTS tree
    Each node own its SMILES, number of visits and score
    """

    def __init__(self, smiles=SMILES(), parent=None):
        """
        Initialisation of a new node

        :param smiles: the state of the SMILES on the node
        :type smiles: SMILES
        :param parent: parent of the node
        :type parent: Node
        """
        self.smiles = smiles
        self.parent = parent
        self.children = []
        self.visits = 0
        self.score = 0.0

    def new_child(self, child_smile):
        """
        Add a new child to the current node

        :param child_smile: the SMILES used to create the new node
        :type child_smile: SMILES
        :return: None
        """
        self.children.append(Node(child_smile, self))

    def update(self, reward):
        """
        Recursive function to update a Node

        :param reward: reward to update the score of the Node
        :type reward: float
        :return: None
        """
        self.score += reward
        self.visits += 1
        if self.parent:
            self.parent.update(reward)

    def fully_expanded(self):
        """
        Return True if the number of children is equal to the size of the vocabulary (minus 2 for \\n and &)
        or if the SMILES of the node is greater than 80 (maximum for the RNN)

        :return: True if the node is terminal in the tree else False
        """
        return (len(self.children) == (len(p.tokens) - 2)) or (len(self.smiles.element) > 80)

    def out_pptree(self, parent=None):
        """
        Recursive function
        Return instance of pptree module to print the tree

        :param parent: parent of the tree (for the recursive part of the function)
        :type parent: Node
        :return: instance of pptree.Node
        """
        name = repr(self)
        if parent:
            current = pptree.Node(name, parent)
        else:
            current = pptree.Node(name)
        for c in self.children:
            c.out_pptree(current)
        return current

    def get_height(self):
        """
        Return the height of the node

        :return: the height of the node
        """
        if not self.children:
            return 1
        else:
            return 1 + max([n.get_height() for n in self.children])

    def get_size(self):
        """
        Return the size of the Node

        :return: the size of the Node
        """
        if not self.children:
            return 1
        else:
            return 1 + sum([n.get_size() for n in self.children])

    def echo(self):
        """
        Print the current tree

        :return: None
        """
        pptree.print_tree(self.out_pptree())
        print("Size of the tree : %d" % self.get_size())
        print("Height of the tree : %d" % self.get_height())

    def __eq__(self, other):
        """
        Equality between two nodes

        :param other: an other node
        :type other: Node
        :return: True if other is the same node
        """
        return str(self.smiles) == str(other.smiles)

    def __repr__(self):
        """
        Representation of a Node

        :return: string representing the node
        """
        return str(self.smiles) + " " + str(int(self.visits)) + " " + str(round(self.score, 2))
