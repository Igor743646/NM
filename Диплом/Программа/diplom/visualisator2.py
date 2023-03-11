import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib import cm 
from scipy.special import gamma as GAMMA, factorial, beta as BETA 
import math

file = open("out.txt", "r")

result1 = []

number_of_strings = file.readline()
h1, tau1, L1 = float(file.readline()), float(file.readline()), float(file.readline())

for i in range(int(number_of_strings)):
    s = file.readline()
    a = s[:-1].strip().split(" ")

    result1.append([float(num) for num in a])

file.close()

decision = lambda x,t : (x**2) * (t**2)

def draw(res, rows, decision = None):

    colums = 4
    all_count = colums * rows
    step = len(res[0]) // (all_count - 1)

    fig, ax = plt.subplots(rows, colums, figsize=(16, 2*rows))

    X_decision = np.linspace(L1, L1 + h1 * len(res[0][0]), 200)
    X1 = np.array([L1 + h1*i for i in range(len(res[0][0]))])
    
    for kk in range(0, all_count):

        TIME_decision = np.array([tau1*kk*step] * len(X_decision))

        if not decision == None:
            ax[kk // colums, kk % colums].plot(X_decision, decision(X_decision, TIME_decision))

        ax[kk // colums, kk % colums].plot(X1, res[0][kk*step], '.-', linewidth = 1)
        ax[kk // colums, kk % colums].set_title(f"{kk*step} (t = {tau1*kk*step})")

    plt.show()


def plot3D(res, decision, vmax = None):

    XX = np.linspace(L1, L1 + h1 * len(res[0]), len(res[0]))
    TT = np.linspace(0.0, tau1 * len(res), len(res))
    XX, TT = np.meshgrid(XX, TT)
    ZZ = decision(XX, TT)
    
    fig = plt.figure(figsize=(16,8))
    ax1 = fig.add_subplot(2, 4, 1, projection = "3d")
    ax2 = fig.add_subplot(2, 4, 2, projection = "3d")
    ax3 = fig.add_subplot(2, 4, (3,4))
    
    ax1.plot_surface(XX, TT, ZZ, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    ax1.view_init(40, -60)
    
    ax2.plot_surface(XX, TT, np.array(res), cmap='viridis', linewidth=0, antialiased=False)
    ax2.view_init(40, -60)
    
    cf = ax3.contourf(XX, TT, np.abs(ZZ - np.array(res)), 200, cmap=cm.coolwarm, vmax = vmax)
    fig.colorbar(cf)
    
    plt.show()

# draw([result1], 6, decision)

plot3D(result1, decision)