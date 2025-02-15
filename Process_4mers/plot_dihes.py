import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Rectangle
import pandas as pd

def plot_scatter(data, file_name, add_guide=False):
    fig, ax = plt.subplots()
    sns.scatterplot(data, x='Phi', y='Psi', size=1, alpha=0.15)
    if add_guide:
        ax.add_patch(Rectangle((180,105), 90, 90, edgecolor='black', fill=False, lw=2))
        plt.annotate(r"$\beta$", (220,140), fontsize=18)
        ax.add_patch(Rectangle((60,-90), 45, 60, edgecolor='black', fill=False, lw=2))
        plt.annotate(r"$\gamma$", (75,-65), fontsize=18)
        ax.add_patch(Rectangle((225,-60), 90, 105, edgecolor='black', fill=False, lw=2))
        plt.annotate(r"$\delta$", (260,-10), fontsize=18)
        ax.add_patch(Rectangle((285,-60), 30, 30, edgecolor='black', fill=False, lw=2))
        plt.annotate(r"$\alpha$", (290,-50), fontsize=18)
        ax.add_patch(Rectangle((45,120), 135, 120, edgecolor='black', fill=False, lw=2))
        plt.annotate(r"$\varepsilon$", (105,180), fontsize=18)
        ax.add_patch(Rectangle((195,45), 60, 60, edgecolor='black', fill=False, lw=2))
        plt.annotate(r"$\zeta$", (220,70), fontsize=18)
        ax.add_patch(Rectangle((270,120), 45, 75, edgecolor='black', fill=False, lw=2))
        plt.annotate(r"$P_{II}$", (285,145), fontsize=18)
        ax.add_patch(Rectangle((255,45), 45, 60, edgecolor='black', fill=False, lw=2))
        plt.annotate(r"$\gamma'$", (270,70), fontsize=18)
        ax.add_patch(Rectangle((30,-15), 90, 90, edgecolor='black', fill=False, lw=2))
        plt.annotate(r"$\delta'$", (70,25), fontsize=18)
        ax.add_patch(Rectangle((45,165), 75, 75, edgecolor='black', fill=False, lw=2))
        plt.annotate(r"$P_{II}'$", (80,200), fontsize=18)
    plt.xlabel('Phi')
    plt.ylabel('Psi')
    l = plt.legend()
    l.remove()
    fig.savefig(file_name)
    plt.close(fig)

df = pd.read_csv('selected_4mers.csv')
phi, psi = [], []

for i, row in df.iterrows():
    for n in [1, 2, 3, 4]:
        phi.append(row[f'Phi {n}'])
        psi.append(row[f'Psi {n}'])

plot_df = pd.DataFrame({'Phi': phi, 'Psi': psi})
plot_scatter(plot_df, 'selected_4mers_dihe.png', True)