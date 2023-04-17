import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib import cm 
from scipy.special import gamma as GAMMA, factorial, beta as BETA 
import math

result1 = []

with open("out.txt", "r") as file:
    number_of_strings = file.readline()
    h1, tau1, L1 = float(file.readline()), float(file.readline()), float(file.readline())

    for i in range(int(number_of_strings)):
        s = file.readline()
        a = s[:-1].strip().split(" ")

        result1.append([float(num) for num in a])


result2 = []

# open("task/canal/canal_map.txt", "r")
with open("task/canal_with_crack/canal_with_crack_map.txt", "r") as file:
    number_of_strings = 3000
    h2, tau2, L2 = 0.005, 1/3000, -0.5

    for i in range(int(number_of_strings)):
        s = file.readline()
        a = s[:-1].strip().split(" ")

        result2.append([float(num) for num in a])

# standart diffusion
def decision(x, t):
    sig2 = lambda t : 4 * 0.02 * t
    res = 0.0047 * np.exp(-(x)**2 / sig2(t)) / np.sqrt(math.pi * sig2(t))
    
    if (type(x[0]) is not np.float64):
        res[0] = np.zeros_like(x[0])
        res[0][len(res[0]) // 2 + 1] = 1.0
    elif (t[0] == 0.0):
        res = np.zeros_like(x)
        res[len(res) // 2] = 1.0
    return res

# anomal diffusion
def decision(x, t):
    sig2 = lambda t : 4 * 0.01 * t
    res = 0.005 * np.exp(-(x)**2 / sig2(t)) / np.sqrt(math.pi * sig2(t))
    
    if (type(x[0]) is not np.float64):
        res[0] = np.zeros_like(x[0])
        res[0][len(res[0]) // 2 + 1] = 1.0
    elif (t[0] == 0.0):
        res = np.zeros_like(x)
        res[len(res) // 2] = 1.0
    return res
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

def plot3D(res, decision, vmax = None):

    XX = np.linspace(L2, L2 + h2 * len(res[0]), len(res[0]))
    TT = np.linspace(0, tau2 * len(res), len(res))
    XX, TT = np.meshgrid(XX, TT)
    ZZ = decision(XX, TT)
    
    fig = plt.figure(figsize=(16,8))
    ax1 = fig.add_subplot(2, 4, (1, 5), projection = "3d")
    ax2 = fig.add_subplot(2, 4, (2, 6), projection = "3d")
    ax3 = fig.add_subplot(2, 4, (3,4))
    ax4 = fig.add_subplot(2, 4, (7,8))

    ax1.plot_surface(XX, TT, ZZ, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    ax1.view_init(40, -60)
    
    ax2.plot_surface(XX, TT, np.array(res), cmap='viridis', linewidth=0, antialiased=False)
    ax2.view_init(40, -60)
    
    cf = ax3.contourf(XX, TT, np.abs(ZZ - np.array(res)), 200, cmap=cm.coolwarm, vmax = vmax)
    fig.colorbar(cf)

    ax4.plot(TT[100:], np.max(np.abs(ZZ - np.array(res))[100:], 1) )

    max_x = np.max(np.abs(ZZ - np.array(res)),1)[100:]
    min_max_x = np.sum(max_x*max_x)
    plt.title(f"eps = {min_max_x}")

def plot3D_num(res, decision, vmax = None):

    XX = np.linspace(L2, L2 + h2 * len(res[0]), len(res[0]))
    TT = np.linspace(0, tau2 * len(res), len(res))
    XX, TT = np.meshgrid(XX, TT)
    
    fig = plt.figure(figsize=(16,8))
    ax1 = fig.add_subplot(2, 4, (1, 5), projection = "3d")
    ax2 = fig.add_subplot(2, 4, (2, 6), projection = "3d")
    ax3 = fig.add_subplot(2, 4, (3,4))
    ax4 = fig.add_subplot(2, 4, (7,8))

    ax1.plot_surface(XX, TT, decision, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    ax1.view_init(40, -60)
    
    ax2.plot_surface(XX, TT, np.array(res), cmap='viridis', linewidth=0, antialiased=False)
    ax2.view_init(40, -60)
    
    cf = ax3.contourf(XX, TT, np.abs(decision - np.array(res)), 200, cmap=cm.coolwarm, vmax = vmax)
    fig.colorbar(cf)

    ax4.plot(TT[100:], np.max(np.abs(decision - np.array(res))[100:], 1) )

    max_x = np.max(np.abs(decision - np.array(res)),1)[100:]
    min_max_x = np.sum(max_x*max_x)
    plt.title(f"eps = {min_max_x}")


draw([result1[:3000], result2], 3, decision)

plot3D(result2, decision)
plot3D_num(result2, np.array(result1[:3000]))

plt.show()