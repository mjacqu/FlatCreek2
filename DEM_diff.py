import rasterio
import matplotlib.pyplot as plt
import numpy as np
import glob
import os

# load ortho as background
datapath = '/Users/mistral/Documents/CUBoulder/Science/Sulzer+/data/ChrisLarsen'
imgs19 = glob.glob(os.path.join(datapath,'Aug30_2019_Sulzer*.tif'))

def make_img_list(fps):
    newlist = []
    for fp in fps:
        src = rasterio.open(fp)
        newlist.append(src)
    return newlist

mosaic_list19 = make_img_list(imgs19)

band_order = [0,1,2]
res = 5
crs = 'EPSG:32607'

#load DEM difference
dem_diff = rasterio.open(
            '/Users/mistral/Documents/CUBoulder/Science/Sulzer+/QGIS_Analysis/2019minus2016_clipped.tif',
            ).read(1)

dem_diff = dem_diff - 0.3 #correct by median value over runout zone
dem_diff[dem_diff < -1000] = np.nan

#x and y limits
glacier_view = []

f, ax = plt.subplots(figsize = (7.5, 9))
ax.imshow(
            dem_diff,
            vmin = -10, vmax = 10,
            cmap = 'seismic_r',
            #extent = [dem_2019_src.bounds.left, dem_2019_src.bounds.right, dem_2019_src.bounds.bottom, dem_2019_src.bounds.top]
)
#ax.set_xlim([473844,474202])
#ax.set_ylim([6839191,6839371])
ax.ticklabel_format(useOffset=None, style = 'plain')
ax.set_xlabel('Northing (m)', fontsize = 16)
ax.set_ylabel('Easting (m)', fontsize = 16)
f.tight_layout()
f.show()

f.clf()
plt.close()
