import numpy as np

class DissSimilarityScore():
    '''
    Kwargs (key word arguments) == scores

    You can use the Kwargs/Scores to input as many scores as you want.

    Each score should have a title and a list/tuple of the sequences to score, so in the form {Head3on5: [Head3, Head5]} corresponding to G(H3|H5)
    If you put in more by accident it will only take the first two.
    '''
    def __init__(self, length: int = 10, print_score = False, *args, **scores) -> None:

        # print(scores, type(scores))
        
        self._weights(length = length)
        # print(self._weight)

        for key, sequences in scores.items():
            # print(f"Key = {key}\nSequences = {sequences}")
            local_g = self._little_g(sequence1 = sequences[0], sequence2 = sequences[1])
            local_score = np.sum(self._weight * local_g)
            # print(local_g)
            # print(local_score)

            setattr(self, key, local_score)

            if print_score:
                print(f"{key}\t-\t{local_g}\t-\t{local_score}")



    def _weights(self, length: int = 10):
        '''
        Weights are caculated as Sum[1/2**i] from i = 1 to L. There is nothing saying you can't do this from 0 instead of 1, but the advantage of using 1 as your start index
        is that the score is scaled from 0 to 1. If you started from 1 the score is scaled from 0 to 2. No change in physical meaning, just adjusts your R value.
        '''
        self._weight = np.array([1/(2**(i + 1)) for i in range(length)])


    def _little_g(self, sequence1: str, sequence2: str):
        '''
        '''
        delta_len = len(sequence1) - len(sequence2)

        if delta_len > 0:
            for _ in range(delta_len):
                sequence2 = f"{sequence2}x"

        if delta_len < 0:
            for _ in range(abs(delta_len)):
                sequence1 = f"{sequence1}x"

        g_vec = np.array([0 if seq1 in sequence2[i] else 1 for i, seq1 in enumerate(sequence1)])

        return g_vec


def main():
    '''
    '''


if __name__ in "__main__":
    main()