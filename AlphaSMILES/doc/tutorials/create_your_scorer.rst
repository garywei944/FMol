Create your scorer
==================

To explore the search tree, the MCTS need a score on each node. this score is calculated with a class *Scorer*.


To create a new scorer you have to code a class extended the *Scorer* class and override the *score* function.

.. code-block:: python

    class ScorerForYourNewScore(Scorer):
    """
    Sub class of Scorer
    Describe your new scorer here
    """

    def __init__(self, alpha=0.5):
        super().__init__(alpha)
        print("Scorer for what you are scoring")

    def score(self, smiles):
        """
        give a score to the SMILES

        :param smiles: SMILES to score
        :return: a number positive or negative
        """
        # some code to find the score based on the calculated properties
        return score

During the execution, this function is call by the reward function from the *Scorer* class.
The reward use the score to reduce it between -1 and 1 with the function :

.. math::

    reward = { | alpha * score | \over 1 + | alpha * score |}

