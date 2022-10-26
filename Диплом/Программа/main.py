import numpy as np
import scipy.linalg
import math
import matplotlib.pyplot as plt

from matplotlib import cm
from scipy.special import gamma as GAMMA, factorial, beta as BETA
from MFDES import *
from FDES import *

class Problem1():

    def __init__(self):
        self.L, self.R = -0.1, 0.1
        self.T = 2.0
        self.alpha = 1.7
        self.gamma = 0.9
        self.h = 0.02
        self.tau = 0.05
        self.beta = 0.0

        self.D = lambda x, t: 0.0005
        self.psi = lambda x: 20.0 * 1.0 if np.abs(x) < self.h else 0.0
        self.f = lambda x, t: 0.0
        self.phiL = None
        self.phiR = None
        self.decision = None

class Problem2():

    def __init__(self):
        self.L, self.R = 0.0, 1.0
        self.T = 2.0
        self.alpha = 1.7
        self.gamma = 0.7
        self.h = 0.05
        self.tau = 0.05
        self.beta = 1.0

        self.D = lambda x, t: GAMMA(3-self.alpha)/GAMMA(3) * (x**self.alpha)
        self.psi = lambda x: 0.0
        self.f = lambda x, t: (GAMMA(3)/GAMMA(3-self.gamma) * (t**(2-self.gamma)) * (x**2)) - (t * x)**2
        self.phiL = lambda t: 0.0
        self.phiR = lambda t: t**2.0
        self.decision = lambda x, t: (t * x)**2

class Problem3():

    def __init__(self):
        self.L, self.R = -0.1, 0.1
        self.T = 2.0
        self.alpha = 1.7
        self.gamma = 0.9
        self.h = 0.002
        self.tau = 0.02
        self.beta = 0.0

        self.D = lambda x, t: 0.0005
        self.psi = lambda x: 1000.0 * 1.0 if np.abs(x) < self.h else 0.0
        self.f = lambda x, t: 0.0
        self.phiL = None
        self.phiR = None
        self.decision = None

if __name__ == "__main__":

    # P = Problem1()

    # a = FDES(P.L, P.R, P.T, P.alpha, P.beta, P.gamma, P.D, P.h, P.tau, P.psi, P.f, P.phiL, P.phiR).solve()
    # a.draw(3, decision = P.decision)
    # a.draw_first_last(decision = P.decision)
    # a.draw3D(decision = P.decision)
    # a.delta(decision = P.decision)

    # a = MFDES(P.L, P.R, P.T, P.alpha, P.beta, P.gamma, P.D, P.h, P.tau, P.psi, P.f, P.phiL, P.phiR).solve()
    # a.draw(3, decision = P.decision)
    # a.draw_first_last(decision = P.decision)
    # a.draw3D(decision = P.decision)
    # a.delta(decision = P.decision)

    # P = Problem2()

    # a = FDES(P.L, P.R, P.T, P.alpha, P.beta, P.gamma, P.D, P.h, P.tau, P.psi, P.f, P.phiL, P.phiR).solve()
    # a.draw(3, decision = P.decision)
    # a.draw_first_last(decision = P.decision)
    # a.draw3D(decision = P.decision)
    # a.delta(decision = P.decision)

    # a = MFDES(P.L, P.R, P.T, P.alpha, P.beta, P.gamma, P.D, P.h, P.tau, P.psi, P.f, P.phiL, P.phiR).solve()
    # a.draw(3, decision = P.decision)
    # a.draw_first_last(decision = P.decision)
    # a.draw3D(decision = P.decision)
    # a.delta(decision = P.decision)

    P = Problem3()

    a = FDES(P.L, P.R, P.T, P.alpha, P.beta, P.gamma, P.D, P.h, P.tau, P.psi, P.f, P.phiL, P.phiR).solve()
    a.draw(3, decision = P.decision)
    a.draw_first_last(decision = P.decision)
    a.draw3D(decision = P.decision)
    a.delta(decision = P.decision)

    a = MFDES(P.L, P.R, P.T, P.alpha, P.beta, P.gamma, P.D, P.h, P.tau, P.psi, P.f, P.phiL, P.phiR).solve()
    a.draw(3, decision = P.decision)
    a.draw_first_last(decision = P.decision)
    a.draw3D(decision = P.decision)
    a.delta(decision = P.decision)