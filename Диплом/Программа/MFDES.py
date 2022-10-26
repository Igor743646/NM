import numpy as np
import scipy.linalg
import math
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.special import gamma as GAMMA, factorial, beta as BETA
from FDESBase import *

# ModifiedFractionalDifferentialEquationSolver
class MFDES(FDESBase):
    
    def __init__(self, L, R, T, alpha, beta, gamma, D, h, tau, psi, f, phiL=None, phiR=None):
        super().__init__(L, R, T, alpha, beta, gamma, D, h, tau, psi, f, phiL, phiR)
        
    def solve(self):
        
        for t in range(1, self.k + 1):
            
            # Заполняем вектор d
            d = np.array([0.0 for i in range(self.n + 1)])
        
            for i in range(self.n + 1):
                for j in range(1, t + 1):
                    d[i] += self.g(self.gamma, j) * self.result[t-j][i]

                d[i] = d[i] - (self.tau**self.gamma) * self.f(self.x_i(i), self.t_k(t))
            
            # Заполняем матрицу А
            for i in range(self.n + 1):
                for j in range(self.n + 1):
                    if i == j:
                        self.A[i][j] = self.a(self.x_i(i), self.t_k(t)) * self.g(self.alpha, 1.0) + \
                                       self.b(self.x_i(i), self.t_k(t)) * self.g(self.alpha, 1.0) - 1.0
                    elif i + 1 < j:
                        self.A[i][j] = self.b(self.x_i(i), self.t_k(t)) * self.g(self.alpha, j - i + 1.0)
                    elif i > j + 1:
                        self.A[i][j] = self.a(self.x_i(i), self.t_k(t)) * self.g(self.alpha, i - j + 1.0)
                    elif i == j + 1:
                        self.A[i][j] = self.a(self.x_i(i), self.t_k(t)) * self.g(self.alpha, 2.0) + \
                                       self.b(self.x_i(i), self.t_k(t))
                    else:
                        self.A[i][j] = self.a(self.x_i(i), self.t_k(t)) + \
                                       self.b(self.x_i(i), self.t_k(t)) * self.g(self.alpha, 2.0)
            
            # Если есть граничные условия
            if not(self.phiL == None) and not(self.phiR == None):
                d[0] = self.phiL(self.t_k(t))
                d[-1] = self.phiR(self.t_k(t))
                
                for i in range(self.n + 1):
                    self.A[0][i] = 1.0 if i == 0 else 0.0
                    self.A[self.n][i] = 1.0 if i == self.n else 0.0
            
            # Решаем систему
            self.result[t] = scipy.linalg.solve(self.A, d)
            
            # Если есть граничные условия
            if not(self.phiL == None) and not(self.phiR == None):
                self.result[t][0] = self.phiL(self.t_k(t))
                self.result[t][-1] = self.phiR(self.t_k(t)) 
            
            
        return self