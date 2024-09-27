import math as mt
import numpy as np
from matplotlib import pyplot as plt
from os import *

# definindo as constantes

C1     = 9.0*(10**(-1))
C2     = 2.0*(10**(-3))
C3     = 2.0*(10**(-2))
C4     = 8.6*(10**(-4))
C5     = 2.3*(10**(-4))
C6     = 8.2*(10**(-3))
C7     = 3.7*(10**(-1))
C8     = 1.0*(10**(0))
C9     = 3.2*(10**(0))
C10    = 3.6*(10**(-1))
C11    = 4.0*(10**(-2))
Bm_max = 5.0*(10**(1))
Bn0    = 5.0*(10**(1))

#definindo as funções
def f_V(V, A):
    return C1*V - C2*V*A

def f_Bn(Bn, V):
    return C3*(Bn0 - Bn) - C4*V*Bn

def f_Ba(Ba, Bm, Bn, V):
    return C4*Bn*V - C5*Ba + C6*Bm*V - C7*Ba

def f_Bm(Bm, Ba, V):
    return C5*Ba + C8*Bm*(1-(Bm/Bm_max)) - C6*Bm*V

def f_A(A, Ba, Bm):
    return C9*Ba + C10*Bm - C11*A

def rkc(V:np.ndarray, Bn:np.ndarray, Ba:np.ndarray, Bm:np.ndarray, A:np.ndarray, dt, t_max):
    n_points = int(t_max/dt)
    
    for i in range(n_points-1):
        # resolvendo vírus
        f1 = f_V(V[i], A[i])
        f2 = f_V(V[i]+(dt/2.0)*f1, A[i])
        f3 = f_V(V[i]+(dt/2.0)*f2, A[i])
        f4 = f_V(V[i]+dt*f3, A[i])
        V[i+1] = V[i] + (dt/6.0)*(f1+2.0*f2+2.0*f3+f4)
        
        # resolvendo célula B naive
        f1 = f_Bn(Bn[i], V[i])
        f2 = f_Bn(Bn[i]+(dt/2.0)*f1, V[i])
        f3 = f_Bn(Bn[i]+(dt/2.0)*f2, V[i])
        f4 = f_Bn(Bn[i]+dt*f3, V[i])
        Bn[i+1] = Bn[i] + (dt/6.0)*(f1+2.0*f2+2.0*f3+f4)
        
        # resolvendo célula B de memória
        f1 = f_Bm(Bm[i], Ba[i], V[i])
        f2 = f_Bm(Bm[i]+(dt/2.0)*f1, Ba[i], V[i])
        f3 = f_Bm(Bm[i]+(dt/2.0)*f2, Ba[i], V[i])
        f4 = f_Bm(Bm[i]+dt*f3, Ba[i], V[i])
        Bm[i+1] = Bm[i] + (dt/6.0)*(f1+2.0*f2+2.0*f3+f4)
        
        # resolvendo célula B de ativa
        f1 = f_Ba(Ba[i], Bm[i], Bn[i], V[i])
        f2 = f_Ba(Ba[i]+(dt/2.0)*f1, Bm[i], Bn[i], V[i])
        f3 = f_Ba(Ba[i]+(dt/2.0)*f2, Bm[i], Bn[i], V[i])
        f4 = f_Ba(Ba[i]+dt*f3, Bm[i], Bn[i], V[i])
        Ba[i+1] = Ba[i] + (dt/6.0)*(f1+2.0*f2+2.0*f3+f4)
        
        # resolvendo Anticorpo
        f1 = f_A(A[i], Ba[i], Bm[i])
        f2 = f_A(A[i]+(dt/2.0)*f1, Ba[i], Bm[i])
        f3 = f_A(A[i]+(dt/2.0)*f2, Ba[i], Bm[i])
        f4 = f_A(A[i]+dt*f3, Ba[i], Bm[i])
        A[i+1] = A[i] + (dt/6.0)*(f1+2.0*f2+2.0*f3+f4)

dt = 0.0001
t_max = 120

n_points = int(t_max/dt)

t = np.arange(0, t_max+(dt/4.0), dt)
t_plot = np.arange(0, 20.0+(dt/4.0), dt)

V  = np.zeros(n_points)
Bn = np.zeros(n_points)
Ba = np.zeros(n_points)
Bm = np.zeros(n_points)
A  = np.zeros(n_points)

V[0]  = 10.0 
Bn[0] = 350.0 
Ba[0] = 0.0 
Bm[0] = 0.0 
A[0]  = 0.0

results_dir = 'python_results'

if not path.isdir(results_dir):
    makedirs(results_dir)

rkc(V=V, Bn=Bn, Ba=Ba, Bm=Bm, A=A, dt=dt, t_max=t_max)

plot_points = n_points/10 

if __name__ == "__main__":
    with plt.style.context('seaborn-v0_8'):

        plt.figure()
        plt.plot(t_plot, V[:len(t_plot)], color='red', linestyle='--', label='Concentração de vírus')
        plt.xlabel("Tempo")
        plt.ylabel("Valor")
        plt.title("Variação do vírus no tempo", fontsize=10)
        plt.legend(loc='best', fontsize=8)
        plt.grid(True)
        plt.tight_layout()  # Minimiza margens e espaços
        plt.savefig(results_dir+"/"+"V.pdf", bbox_inches='tight', pad_inches=0.05)
        plt.close()

        plt.figure()
        plt.plot(t_plot, Bm[:len(t_plot)], color='green', linestyle='--', label='Concentração de Bm')
        plt.xlabel("Tempo")
        plt.ylabel("Valor")
        plt.title("Variação de B de memória no tempo", fontsize=10)
        plt.legend(loc='best', fontsize=8)
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(results_dir+"/"+"Bm.pdf", bbox_inches='tight', pad_inches=0.05)
        plt.close()

        plt.figure()
        plt.plot(t_plot, Bn[:len(t_plot)], color='magenta', linestyle='--', label='Concentração de Bn')
        plt.xlabel("Tempo")
        plt.ylabel("Valor")
        plt.title("Variação de B de naive no tempo", fontsize=10)
        plt.legend(loc='best', fontsize=8)
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(results_dir+"/"+"Bn.pdf", bbox_inches='tight', pad_inches=0.05)
        plt.close()

        plt.figure()
        plt.plot(t_plot, Ba[:len(t_plot)], color='purple', linestyle='--', label='Concentração de Ba')
        plt.xlabel("Tempo")
        plt.ylabel("Valor")
        plt.title("Variação de B de ativa no tempo", fontsize=10)
        plt.legend(loc='best', fontsize=8)
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(results_dir+"/"+"Ba.pdf", bbox_inches='tight', pad_inches=0.05)
        plt.close()

        plt.figure()
        plt.plot(t[:len(t)-1], A, color='orange', linestyle='--', label='Concentração de A')
        plt.xlabel("Tempo")
        plt.ylabel("Valor")
        plt.title("Variação de Anticorpo no tempo", fontsize=10)
        plt.legend(loc='best', fontsize=8)
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(results_dir+"/"+"A.pdf", bbox_inches='tight', pad_inches=0.05)
        plt.close()
