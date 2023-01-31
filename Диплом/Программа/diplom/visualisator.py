import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib import cm 
from scipy.special import gamma as GAMMA, factorial, beta as BETA 

file = open("out.txt", "r")

result = []

for i in file:
    a = i[:-1].strip().split(" ")

    result.append([float(i) for i in a])

file.close()

h = tau = 0.05
L = 0.0

decision = lambda x, t: (t * x)**2.0

def draw(rows, decision = None):
    res = result
    colums = 3
    all_count = colums * rows
    step = len(result) // (all_count - 1)

    fig, ax = plt.subplots(rows, colums, figsize=(8, 2*rows))
    X = np.array([L + h*i for i in range(len(result[0]))])

    for kk in range(0, all_count):

        TIME = np.array([tau*kk*step for i in range(len(result[0]))])
        
        if not decision == None:
            ax[kk // colums, kk % colums].plot(X, decision(X, TIME))
        ax[kk // colums, kk % colums].plot(X, res[kk*step], '.-')
        ax[kk // colums, kk % colums].set_title(f"{kk*step} (t = {tau*kk*step})")

    plt.show()

draw(3, decision)