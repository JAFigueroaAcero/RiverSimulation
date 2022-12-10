import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator

def e(c,amp,x):
    return np.exp(-(x-c)**2/(2*amp**2))/ (amp * np.sqrt(2 * np.pi))

def S(x,t,i):
    b = i*(99461*e(3470,4399,x) + 56729 * e(21730,3500,x) + 1190 * e(33990,600,x) + 400 * e(53690,200,x))
    return (1/(24))*(1-np.cos(np.pi/(12) * t))*b


def main():
    Ttotal = 80
    dt = 0.05
    dx = 1000
    L = 200000
    mu = 0.26
    u = 0.5*3600
    q = 0.01
    i = 0.131
    xlen = int(L/dx)
    tlen = int(Ttotal/dt)
    C = [[0 for a in range(xlen)]]
    for t in range(1,tlen):
        loc_l = []
        loc_1 = dt * (2*mu * (C[t-1][1]-C[t-1][0]- 2*dx*q)/(dx**2) - u*q + S(0,t*dt,i)) + C[t-1][0]
        if loc_1 < 0:
            loc_l.append(0)
        else:
            loc_l.append(loc_1)
        
        for x in range(1,xlen-1):
            loc_2 = dt * (mu * (C[t-1][x+1]-2*C[t-1][x] + C[t-1][x-1])/(dx**2) - u*(C[t-1][x+1]-C[t-1][x-1])/(2*dx) + S(x*dx,t*dt,i)) + C[t-1][x]
            if loc_2 < 0:
                loc_l.append(0)
            else:
                loc_l.append(loc_2)
        loc_l.append(0)
        C.append(loc_l)
    C = np.array([c[0:60] for c in C])
    X = np.arange(0,60000/1000,dx/1000)
    Y = np.arange(0,Ttotal,dt)
    X,Y = np.meshgrid(X,Y)
    print(X.shape,Y.shape,C.shape)
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    surf = ax.plot_surface(X, Y, C,cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter('{x:.02f}')
    ax.set_xlabel(r'$Distancia-(km)$')
    ax.set_ylabel(r'$Tiempo-(h)$')
    ax.set_zlabel(r'$Contaminante - \left(\frac{kg}{L}\right)$')
    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.show()
    

if __name__ == '__main__':
    main()