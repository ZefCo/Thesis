import numpy as np


class DistanceMatrix():
    def __init__(self, x: np.array, method = "euclidean"):
        self.x: np.ndarray = x
        self.l = self.x.shape[0]
        self.p = int(self.points(self.l))

        if method in "manhattan":
            self.func = self.manhattan
        elif method in "esquare":
            self.func = self.euclidean_square
        else:
            self.func = self.euclidean
        
        self.compute_dist()


    def manhattan(self, xy1, xy2):
        '''
        '''
        return abs(xy1[0] - xy2[0]) + abs(xy1[1] - xy2[1])


    def euclidean_square(self, xy1, xy2):
        '''
        '''
        return (xy1[0] - xy2[0])**2 + (xy1[1] - xy2[1])**2


    def euclidean(self, xy1, xy2):
        '''
        '''
        return np.sqrt((xy1[0] - xy2[0])**2 + (xy1[1] - xy2[1])**2)


    def points(self, l: float):
        '''
        '''
        return (0.5*(l**2)) - (0.5*l)


    def _init_zero(self, d = 1):
        '''
        '''

        try:
            if d == 1:
                self.D1 = np.zeros(self.p)
            elif d == 2:
                self.D2 = np.zeros(shape = (self.l, self.l))

        except np.core._exceptions._ArrayMemoryError as e:
            print("Not enough memory")
            print(f"\tSize = {self.p}\tItem Size = 8-16\tGB = {self.p*8*1E-9} - {self.p*16*1E-9}")
            
            self.D1 = None

        except Exception as e:
            print("New Error when making the zeros matrix")
            print(type(e))
            print(e)
            print(f"\tSize = {self.p}\tItem Size = 8-16\tGB = {self.p*8*1E-9} - {self.p*16*1E-9}")
            
            self.D1 = None


    def compute_dist(self):
        '''
        '''
        self._init_zero()

        if isinstance(self.D1, np.ndarray):
            i = 0
            for k in range(self.l - 1):
                for kk in range(k + 1, self.l):
                    d = self.func(self.x[k], self.x[kk])
                    self.D1[i] = d
                    i += 1

            return self.D1
        
        else:
            return None
        

    def matrix2d(self):
        '''
        '''
        self._init_zero(d = 2)

        for k in range(self.l):
            for i in range(k + 1, self.l):
                self.D2[k][i] = self.D2[i][k] = self.D1[k]

        return self.D2



    # def dis_matrix(matrix) -> np.array:
    #     '''
    #     '''

    #     try:
    #         D = np.zeros(int(p))

    #     except np.core._exceptions._ArrayMemoryError as e:
    #         print("\tStill have a memory error")
    #         print(f"\tSize = {p}\tItem Size = 8?\tGB = {p*8*1E-9}")
    #         return None

    #     except Exception as e:
    #         print(type(e))
    #         print(e)

    #         print(f"\tSize = {p}\tItem Size = 8?\tGB = {p*8*1E-9}")
    #         return None

    #     else:
    #         print(f"\tSize = {D.size}\tItem Size = {D.itemsize}\tGB = {D.size*D.itemsize*1E-9}")

    #         w = len(x)
    #         # print(w)
    #         over_index = 0
    #         for k in range(w):
    #             step = 1 + k
    #             for i in range(0, l - k - 1):
    #                 d = euc(x[i], x[i + step])
    #                 if d < 1:
    #                     print(f"D{i}{i+step}: {d}")


    #                 D[over_index] = d
    #                 over_index += 1

    #         # print(k, w)

    #         return D

