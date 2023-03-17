import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib import cm 
from scipy.special import gamma as GAMMA, factorial, beta as BETA 
import math

file = open("out.txt", "r")

result1 = []
result2 = []


number_of_strings = file.readline()
h1, tau1, L1 = float(file.readline()), float(file.readline()), float(file.readline())

for i in range(int(number_of_strings)):
    s = file.readline()
    a = s[:-1].strip().split(" ")

    result1.append([float(num) for num in a])

number_of_strings = file.readline()
h2, tau2, L2 = float(file.readline()), float(file.readline()), float(file.readline())

for i in range(int(number_of_strings)):
    s = file.readline()
    a = s[:-1].strip().split(" ")

    result2.append([float(num) for num in a])

file.close()

# h, tau = [0.001, 0.001], 0.001

# L = -0.1

sig2 = lambda t : 4 * 0.0005 * t

def decision(x, t):
    return 0.1 * np.exp(-(x)**2 / sig2(t)) / np.sqrt(math.pi * sig2(t)) if t[0] > 0.0 else (10.0 * (np.abs(x) < np.array([h1] * len(x))) )

# decision = lambda x,t : (np.sin(x)) * np.exp(2*t)
# decision = lambda x,t : (x * x) * t * t
# decision = lambda x,t : (x * x) * np.exp(2*t)
# print(decision(0, tau2))
# decision = None

def draw(res, rows, decision = None):
    # if res is not list:
    #     res = [res]

    colums = 4
    all_count = colums * rows
    step = len(res[0]) // (all_count - 1)

    fig, ax = plt.subplots(rows, colums, figsize=(16, 2*rows))

    X_decision = np.linspace(L1, L1 + h1 * len(res[0][0]), 200)
    X1 = np.array([L1 + h1*i for i in range(len(res[0][0]))])
    X2 = np.array([L2 + h2*i for i in range(len(res[1][0]))])
    
    for kk in range(0, all_count):

        # TIME1 = np.array([tau1*kk*step for i in range(len(res[0][0]))])
        # TIME2 = np.array([tau2*kk*step for i in range(len(res[1][0]))])
        TIME_decision = np.array([tau1*kk*step] * len(X_decision))

        if not decision == None:
            ax[kk // colums, kk % colums].plot(X_decision, decision(X_decision, TIME_decision))
        # for solution in range(len(res)):
        ax[kk // colums, kk % colums].plot(X1, res[0][kk*step], '.-', linewidth = 1)
        ax[kk // colums, kk % colums].plot(X2, res[1][kk*step], '.-', linewidth = 1)
        ax[kk // colums, kk % colums].set_title(f"{kk*step} (t = {tau1*kk*step})")

# draw([result1, result2], int(number_of_strings) // 4, decision)
draw([result1, result2], 5, decision)

plt.show()