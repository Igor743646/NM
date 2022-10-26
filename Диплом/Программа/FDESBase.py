import numpy as np
import scipy.linalg
import math
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.special import gamma as GAMMA, factorial, beta as BETA

class FDESBase():

    def __init__(self, L, R, T, alpha, beta, gamma, D, h, tau, psi, f, phiL=None, phiR=None):
        self.n = int ((R - L) / h)
        self.k = int (T / tau)
        print(self.n, self.k)
        self.alpha = alpha
        self.gamma = gamma
        self.tau = tau
        self.h = h
        self.L = L
        self.D = D
        self.f = f
        self.phiL = phiL
        self.phiR = phiR
        self.result = np.array([[0.0 for i in range(self.n + 1)] for t in range(self.k + 1)])
        
        for i in range(self.n + 1):
            self.result[0][i] = psi(self.x_i(i))
            
        self.A = np.array([[0.0 for i in range(self.n + 1)] for t in range(self.n + 1)], dtype=np.float64)
        
        self.a = lambda x, t: (1.0 + beta) * (self.D(x, t) / 2.0) * ((tau**gamma)/(h**alpha))
        self.b = lambda x, t: (1.0 - beta) * (self.D(x, t) / 2.0) * ((tau**gamma)/(h**alpha))

    def g(self, a, j):
        return (-1)**j / (a + 1) / BETA(j + 1, a - j +1)
        
    def x_i(self, i):
        return self.L + i * self.h
    
    def t_k(self, k):
        return k * self.tau

    def draw(self, rows, decision = None):
        res = self.result
        colums = 3
        all_count = colums * rows
        step = self.k // (all_count - 1)

        fig, ax = plt.subplots(rows, colums, figsize=(8, 2*rows))
        X = np.array([self.L + self.h*i for i in range(self.n+1)])

        for kk in range(0, all_count):

            TIME = np.array([self.tau*kk*step for i in range(self.n+1)])
            
            if not decision == None:
                ax[kk // colums, kk % colums].plot(X, decision(X, TIME))
            ax[kk // colums, kk % colums].plot(X, res[kk*step], '.-')
            ax[kk // colums, kk % colums].set_title(f"{kk*step} (t = {self.tau*kk*step})")

        plt.show()
        
    def delta(self, decision = None):
        if decision == None: 
            return

        res = np.array(self.result)
        time, delt = [], []

        _, _ = plt.subplots(1, 1, figsize=(8, 4))

        X = np.array([self.h*_i for _i in range(self.n+1)])

        for i in range(self.k):
            TIME = np.array([self.tau*i for _i in range(self.n+1)])
            aaa = np.abs(res[i] - decision(X, TIME))

            time.append(self.tau*i)
            delt.append(np.max(aaa))

        plt.plot(time, delt)
        plt.show()
        
    def draw_first_last(self, decision = None):
        res = self.result
        
        fig, ax = plt.subplots(1, 1, figsize=(8, 4))
        X = np.array([self.L + self.h*i for i in range(self.n+1)])

        TIME1 = np.array([self.tau*0 for i in range(self.n+1)])
        TIME2 = np.array([self.tau*self.k for i in range(self.n+1)])
        if not decision == None:
            ax.plot(X, decision(X, TIME1))
            ax.plot(X, decision(X, TIME2))
        ax.plot(X, res[0], '.-')
        ax.plot(X, res[-1], '.-')

        plt.show()
        
    def draw3D(self, decision = None):
        XX = np.array([self.x_i(i) for i in range(self.n + 1)])
        TT = np.array([self.t_k(j) for j in range(self.k + 1)])
        XX, TT = np.meshgrid(XX, TT)
        
        if decision != None:
            ZZ = decision(XX, TT)
            fig = plt.figure(figsize=(10,4))
            ax1 = fig.add_subplot(1, 4, 1, projection = "3d")
            ax2 = fig.add_subplot(1, 4, 2, projection = "3d")
            ax3 = fig.add_subplot(1, 4, (3,4))
            ax1.plot_surface(XX, TT, ZZ, cmap=cm.coolwarm, linewidth=0, antialiased=False)
            ax1.view_init(40, +30)

            ax2.plot_surface(XX, TT, np.array(self.result), cmap='viridis', linewidth=0, antialiased=False)
            ax2.view_init(40, +30)

            cf = ax3.contourf(XX, TT, np.abs(ZZ - np.array(self.result)), 200, cmap=cm.coolwarm)
            fig.colorbar(cf)
        else:
            fig = plt.figure(figsize=(8,4))
            ax1 = fig.add_subplot(1, 4, (1,2), projection = "3d")
            ax2 = fig.add_subplot(1, 4, (3,4), projection = "3d")
            
            ax1.plot_surface(XX, TT, np.array(self.result), cmap='viridis', linewidth=0, antialiased=False)
            ax1.view_init(40, -30)
            
            ax2.plot_surface(XX, TT, np.array(self.result), cmap='viridis', linewidth=0, antialiased=False)
            ax2.view_init(40, +30)
        
        plt.show()