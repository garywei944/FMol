import math
from abc import ABC, abstractmethod

from mcts import parameters as p


class Scorer(ABC):
    """
    Abstract class for the calculation of the score of a SMILES
    """

    def __init__(self, alpha=0.5):
        """
        Initialise the Scorer

        :param alpha: number used in the reward method
        """
        super().__init__()
        self.alpha = alpha

    @abstractmethod
    def score(self, smiles):
        """
        Abstact method to calculate the score

        :param smiles: the SMILES object to score
        :return: the score of the SMILES
        """
        pass

    def reward(self, smiles, already):
        """
        Method giving the reward depending on the score of the SMILES and if it have already been calculated

        :param smiles: the SMILES to reward
        :param already: True if the SMILES score have already been calculated
        :return: the reward of the SMILES
        """
        score = self.score(smiles)
        reward = ((self.alpha * score) / (1 + math.fabs(self.alpha * score)))
        return reward if not already or not smiles[p.s_valid] else (reward / 2)


class ScorerValidSMILES(Scorer):
    """
    Sub class of Scorer
    Used to reward on the validity of the generated SMILES
    """

    def __init__(self, alpha=0.5):
        super().__init__(alpha)
        print("Scorer for valid SMILES")

    def score(self, smiles):
        """
        Score the SMILES

        :param smiles: SMILES to score
        :return: 1 if the SMILES if valid else -1
        """
        return 1 if smiles[p.s_valid] else -1


class ScorerDFT(Scorer):
    """
    Sub class of scorer
    Used to reward the DFT score of a SMILES
    """

    def __init__(self, alpha=0.5):
        super().__init__(alpha)
        print("Scorer for DFT")

    def score(self, smiles):
        """
        Score the SMILES

        :param smiles: SMILES to score
        :type smiles: SMILES
        :return: -100 if the SMILES is not valid, -0.5 if the DFT calcul didn't end a value above 0 depending on the peak of oscillator strength otherwise
        """
        if smiles[p.s_valid]:
            try:
                dft = smiles[p.s_dft]
            except Exception as e:
                print(e)
                return -0.5
            score = 0
            for line in dft:
                if 300 <= line["nm"] < 400:
                    score += line["f"] * 25
                if 400 <= line["nm"] < 500:
                    score += line["f"] * 50
                if 500 <= line["nm"]:
                    score += line["f"] * 100
            return score
        else:
            return -100

