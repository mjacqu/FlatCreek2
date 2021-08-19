import os
import pandas as pd
import re
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import ticker
import numpy as np


def read_from_xyz(path, filename):
    with open(os.path.join(path, filename), "r") as file:
        cnt = 0
        pos = 0
        for line in file:
            cnt += 1
            if re.match("/        X       Elevation    Resistivity  Conductivity", line):
                pos = cnt
    data = pd.read_csv(os.path.join(path, filename), header = None, comment = '/',
    sep = '\s+', skiprows = pos, names = ['X', 'Elevation', 'Resistivity', 'Conductivity'])
    return data


def plot_ert(data,filename, ncolours, cmap, clevels = None, save = False):
    if clevels is None:
        clevels = np.logspace(np.log2(np.min(data['Resistivity'])),np.log2(np.max(data['Resistivity'])),num=ncolours,base=2)
    else:
        clevels = np.logspace(np.log2(clevels[0]),np.log2(clevels[1]),num=ncolours,base=2)
    fig, ax = plt.subplots()
    triang = mpl.tri.Triangulation(data['X'], data['Elevation'])
    mask = mpl.tri.TriAnalyzer(triang).get_flat_tri_mask(.1)
    triang.set_mask(mask)
    cc = ax.tricontourf(triang,data['Resistivity'],levels=clevels, norm=mpl.colors.LogNorm(), cmap=cmap)
    #ax.axis('equal') # make axes equal (i.e. 1m x = 1m y)
    cb = fig.colorbar(cc, format='%.0f', orientation = 'horizontal')
    cb.ax.minorticks_off()
    cb.set_label('Resistivity ($\Omega\cdot$m)')
    ax.set_ylim([np.floor(data['Elevation'].min()),np.ceil(data['Elevation'].max())+2])
    fig.gca().set_aspect('equal', adjustable='box')
    ax.set_title(filename)
    ax.set_xlabel('Distance (m)')
    ax.set_ylabel('Depth (m)')
    fig.tight_layout()
    if save is True:
        fig.savefig(filename[:-4] + ".pdf")
    else:
        fig.show()

path = '/Users/mistral/Documents/CUBoulder/Science/Sulzer+/data/ERT_data_from_Matthias/ERTxyz'
filenames = os.listdir(path)
filenames_deposit = ['flatcreek_line_2.xyz', 'flatcreek_line_3.xyz', 'flatcreek_line_4.xyz', 'flatcreek_line_5.xyz',
    'flatcreek_line_6.xyz', 'flatcreek_line_7.xyz', 'flatcreek_line_8.xyz', 'flatcreek_line_9.xyz',
    'eastdeposit_line_11.xyz', 'eastdeposit_line_12.xyz','eastdeposit_line_13.xyz']
deposit_clevels = [30, 15000]
filenames_ice = ['eastdeposit_line_10.xyz']
ice_clevels = [100, 100000]

ncolours=36
cmap = 'magma'

for f in filenames:
    if not f.startswith('.'):
        data = read_from_xyz(path, f)
        plot_ert(data, f, ncolours, cmap, save = True)

for f in filenames_deposit:
    if not f.startswith('.'):
        data = read_from_xyz(path, f)
        plot_ert(data, f, ncolours, cmap, clevels = deposit_clevels, save = True)

for f in filenames_ice:
    if not f.startswith('.'):
        data = read_from_xyz(path, f)
        plot_ert(data, f, ncolours, cmap, clevels = ice_clevels, save = True)
