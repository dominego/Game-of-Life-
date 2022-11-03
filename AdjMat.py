import numpy as np
import random


class AdjMat:
    def __init__(
        self,
        l=6,
        conf=[],
        Mat=[],
        result=1,
    ):
        def SquaredGridAdjMat(self):
            Mat = np.zeros((self.l * self.l, self.l * self.l))
            for i in range(self.l * l):
                Mat[i, (i + 1) % self.l + self.l * (i // self.l)] = 1
                Mat[i, (i + self.l) % (self.l ** 2)] = 1
                self.Mat += np.transpose(Mat)

        def TriangularGridAdjMat(self):
            Mat = np.zeros((self.l * self.l, self.l * self.l))
            for i in range(self.l * l):
                Mat[i, (i + 1) % self.l + self.l * (i // self.l)] = 1
                Mat[i, (i + self.l) % (self.l ** 2)] = 1
            for i in range(l):
                for j in range(0, l, 2):
                    if i != self.l - 1:
                        Mat[(j * self.l) + i, (j * self.l) + i - self.l + 1] = 1
                        Mat[(j * self.l) + i, (j * self.l) + i + self.l + 1] = 1
                    elif i == l - 1:
                        Mat[(j * self.l) + i, (j * self.l) + i + 1 - 2 * self.l] = 1
                        Mat[(j * self.l) + i, (j * self.l) + i + 1] = 1

            self.Mat += np.transpose(Mat)

        def InitializeRandomConfiguration(self):
            for i in range(self.l * self.l):
                self.conf.append(random.randint(0, 1))
        
        def ArbitraryFunction(self):
            if self.Mat!=[] and self.conf!=[]:
                for i in range(self.l * self.l):
                    for j in range(self.l * self.l):
                        if self.Mat[i, j] != 0:
                            self.result *= function(self.conf[i], self.conf[j])
            else:
                return -1

