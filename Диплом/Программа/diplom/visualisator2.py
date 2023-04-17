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

decision = lambda x,t : (x**2) * (t**2)



sig2 = lambda t : 4 * 0.0005 * t

def decision(x, t):
    res = 0.1 * np.exp(-(x)**2 / sig2(t)) / np.sqrt(math.pi * sig2(t)) 
    # res[0] = 10.0 * (np.abs(x[0]) < np.array([h1] * len(x[0])))
    return res

decision = lambda x,t : (np.sin(x)) * np.exp(2*t)
# decision = lambda x,t : (x * x) * np.exp(2*t)
# decision = lambda x,t : (x * x) * t * t
# decision = lambda x,t : t * t
# decision = lambda x, t: np.zeros_like(x)
# decision = None

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



def plot3D(res, decision, vmax = None):

    XX = np.linspace(L1, L1 + h1 * len(res[0]), len(res[0]))
    TT = np.linspace(0, tau1 * len(res), len(res))
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

    ax4.plot(TT, np.max(np.abs(ZZ - np.array(res)), 1) )

def delta12(res1, res2, f):
    XX = np.linspace(L1, L1 + h1 * len(res1[0]), len(res1[0]))
    TT = np.linspace(0, tau1 * len(res1), len(res1))
    XX, TT = np.meshgrid(XX, TT)
    
    fig = plt.figure(figsize=(16,8))
    ax1 = fig.add_subplot(1, 1, 1, projection = "3d")

    ax1.plot_surface(XX, TT, np.array(res1) - np.array(res2) + f(XX, TT) * np.power(TT, 0.8), cmap=cm.coolwarm, linewidth=0, antialiased=False)
    # ax1.plot_surface(XX, TT, f(XX, TT))
    ax1.view_init(40, -60)
    

# draw([result1], 6, decision)

plot3D(result1, decision)
plot3D(result2, decision)
# delta12(result1, result2, f)

plt.show()