import sys
sys.path.append('../ISCE')
import insarhelpers
import geopandas as gpd
import matplotlib.pyplot as plt
import glob
import os
import rasterio.plot
import numpy as np

#Make plot with three tiles showing
file_paths= ('/Users/mistral/Documents/CUBoulder/Science/Sulzer+/data/planet/20120815/files',
            '/Users/mistral/Documents/CUBoulder/Science/Sulzer+/data/planet/20130811/files',
            '/Users/mistral/Documents/CUBoulder/Science/Sulzer+/data/planet/20150813/files')


imgs = [glob.glob(os.path.join(p,'*_SR_clip.tif')) for p in file_paths]

mosaic_list = []
for l in imgs:
    m = []
    for fp in l:
        src = rasterio.open(fp)
        m.append(src)
    mosaic_list.append(m)

band_order = [3,1,0]

pre, extent = insarhelpers.make_MS_img(mosaic_list[0], band_order, brightness_f = 4)
post13, extent = insarhelpers.make_MS_img(mosaic_list[1], band_order, brightness_f = 3)
post15, extent = insarhelpers.make_MS_img(mosaic_list[2], band_order, brightness_f = 4)

#import shapefiles
runout2013 = gpd.read_file('/Users/mistral/Documents/CUBoulder/Science/Sulzer+/QGIS_Analysis/2013_runout.geojson')
runout2015 = gpd.read_file('/Users/mistral/Documents/CUBoulder/Science/Sulzer+/QGIS_Analysis/2015_runout.geojson')
flowlline_2013 = gpd.read_file('/Users/mistral/Documents/CUBoulder/Science/Sulzer+/QGIS_Analysis/2013_runout_distance.geojson')
flowlline_2015 = gpd.read_file('/Users/mistral/Documents/CUBoulder/Science/Sulzer+/QGIS_Analysis/2015_runout_distance.geojson')


x_limits = [469000,479500]
y_limits = [6833500,6844000]

f, axs = plt.subplots(1,2, figsize = (12,6))
#axs[0].imshow(pre, extent = extent)
#axs[0].set_xlim(x_limits)
#axs[0].set_ylim(y_limits)
axs[0].imshow(post13, extent = extent)
runout2013.exterior.plot(ax = axs[0], color = 'khaki', linewidth = 0.3)
flowlline_2013.plot(ax = axs[0], linestyle= '--', color = 'w', linewidth = 0.3)
axs[0].set_xlim(x_limits)
axs[0].set_ylim(y_limits)
axs[1].imshow(post15, extent = extent)
runout2015.exterior.plot(ax = axs[1], color = 'khaki', linewidth = 0.3)
flowlline_2015.plot(ax = axs[1],linestyle= '--', color = 'w', linewidth = 0.3)
axs[1].set_xlim(x_limits)
axs[1].set_ylim(y_limits)
f.tight_layout()
#f.show()
#f.savefig('FlatCreekRunouts.pdf')
