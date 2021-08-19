from windrose import WindroseAxes
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np

data = pd.read_csv('/Users/mistral/Documents/CUBoulder/Science/Sulzer+/data/SedimentAnalysis/sed_orientation.csv')

site1_dir = []
site2_dir = []
site3_dir = []

for i in range(0,len(data)-2):
    tmp = np.ones(data.site_1[i])*data.site_1_2_deg[i]
    site1_dir = np.concatenate([site1_dir,tmp])
    site1_num =np.ones(len(site1_dir))
    tmp = np.ones(data.site_2[i])*data.site_1_2_deg[i]
    site2_dir = np.concatenate([site2_dir,tmp])
    site2_num =np.ones(len(site2_dir))
    tmp = np.ones(data.site_3[i])*data.site_3_deg[i]
    site3_dir = np.concatenate([site3_dir,tmp])
    site3_num =np.ones(len(site3_dir))

dir = site1_dir
num = site1_num
no_dim = data.site_1[12]
vert = data.site_1[13]

ax = WindroseAxes.from_ax()
vert_circ = plt.Circle((0, 0), vert, transform=ax.transData._b, color='r', fill = False, linewidth = 2)
ax.add_artist(vert_circ)
ax.bar(dir, num, bins = 1, opening = 0.95, nsector = 12, color = 'black')
no_dim_circ = plt.Circle((0, 0), no_dim, transform=ax.transData._b, color='blue', fill = False, linewidth = 2)
ax.add_artist(no_dim_circ)
ax.set_yticks(np.arange(10, 60, step=10))
ax.set_yticklabels(np.arange(10, 60, step=10))
plt.show()
#plt.savefig('Orientations_site_1.pdf', bbox_inches = 'tight')
