import numpy as np
from matplotlib import pyplot as plt

def read_file(file_name):
    with open(file_name, 'r') as file_in:
        lines = file_in.readlines()
    array = np.array([float(line.strip()) for line in lines])
    return array

file_name = 'cpp_results/V.txt'
V = read_file(file_name)

file_name = 'cpp_results/Bn.txt'
Bn = read_file(file_name)

file_name = 'cpp_results/Ba.txt'
Ba = read_file(file_name)

file_name = 'cpp_results/Bm.txt'
Bm = read_file(file_name)

file_name = 'cpp_results/A.txt'
A = read_file(file_name)

t = np.linspace(0, 120, len(A))
dt = t[1] - t[0]

t_plot = np.arange(0, 20.0 + (dt / 4.0), dt)

results_dir = 'cpp_results'

print("read...")

if __name__ == "__main__":
    with plt.style.context('seaborn-v0_8-dark'):
        plt.figure()
        plt.plot(t_plot, V[:len(t_plot)], color='red', linestyle='-', label='Concentração de vírus')
        plt.xlabel("Tempo")
        plt.ylabel("Valor")
        plt.title("Variação do vírus no tempo", fontsize=10)
        plt.legend(loc='best', fontsize=8)
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(results_dir + "/V.pdf", bbox_inches='tight', pad_inches=0.05)
        plt.close()

        plt.figure()
        plt.plot(t_plot, Bm[:len(t_plot)], color='green', linestyle='-', label='Concentração de Bm')
        plt.xlabel("Tempo")
        plt.ylabel("Valor")
        plt.title("Variação de B de memória no tempo", fontsize=10)
        plt.legend(loc='best', fontsize=8)
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(results_dir + "/Bm.pdf", bbox_inches='tight', pad_inches=0.05)
        plt.close()

        plt.figure()
        plt.plot(t_plot, Bn[:len(t_plot)], color='magenta', linestyle='-', label='Concentração de Bn')
        plt.xlabel("Tempo")
        plt.ylabel("Valor")
        plt.title("Variação de B de naive no tempo", fontsize=10)
        plt.legend(loc='best', fontsize=8)
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(results_dir + "/Bn.pdf", bbox_inches='tight', pad_inches=0.05)
        plt.close()

        plt.figure()
        plt.plot(t_plot, Ba[:len(t_plot)], color='purple', linestyle='-', label='Concentração de Ba')
        plt.xlabel("Tempo")
        plt.ylabel("Valor")
        plt.title("Variação de B de ativa no tempo", fontsize=10)
        plt.legend(loc='best', fontsize=8)
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(results_dir + "/Ba.pdf", bbox_inches='tight', pad_inches=0.05)
        plt.close()

        plt.figure()
        plt.plot(t, A, color='orange', linestyle='-', label='Concentração de A')
        plt.xlabel("Tempo")
        plt.ylabel("Valor")
        plt.title("Variação de Anticorpo no tempo", fontsize=10)
        plt.legend(loc='best', fontsize=8)
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(results_dir + "/A.pdf", bbox_inches='tight', pad_inches=0.05)
        plt.close()
