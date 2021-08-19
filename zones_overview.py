import numpy as np
import rasterio
import earthpy
import earthpy.plot
import earthpy.spatial
import matplotlib.pyplot as plt
import geopandas as gpd

dem = '/Users/mistral/Documents/CUBoulder/Science/Sulzer+/data/DEM/regionalIFSAR_UTM_10m.tif'
with rasterio.open(dem) as src:
    elevation = src.read(1)
    elevation[elevation < 0] = np.nan

dem_extent = [373271.726, 521821.726, 6743453.939, 6894553.939]
hillshade = earthpy.spatial.hillshade(elevation, azimuth = 135, altitude = 45)

#import outlines:
outline2015 = gpd.read_file('/Users/mistral/Documents/CUBoulder/Science/Sulzer+/QGIS_Analysis/2015_runout.geojson')
hummock_area = gpd.read_file('/Users/mistral/Documents/CUBoulder/Science/Sulzer+/QGIS_Analysis/hummocks_broad.geojson')
levee2013 = gpd.read_file('/Users/mistral/Documents/CUBoulder/Science/Sulzer+/QGIS_Analysis/2013_levee.geojson')
levees = gpd.read_file('/Users/mistral/Documents/CUBoulder/Science/Sulzer+/QGIS_Analysis/levees.geojson')
flowbands = gpd.read_file('/Users/mistral/Documents/CUBoulder/Science/Sulzer+/QGIS_Analysis/pressure_ridges.geojson')
river2016 = gpd.read_file('/Users/mistral/Documents/CUBoulder/Science/Sulzer+/QGIS_Analysis/river.geojson')
river2019= gpd.read_file('/Users/mistral/Documents/CUBoulder/Science/Sulzer+/QGIS_Analysis/river2019.geojson')
detachment16 = gpd.read_file('/Users/mistral/Documents/CUBoulder/Science/Sulzer+/QGIS_Analysis/2016_detachment.geojson')
detachment15 = gpd.read_file('/Users/mistral/Documents/CUBoulder/Science/Sulzer/VolumeCalculations/MainLoss2015.shp')
detachment15 = detachment15.to_crs('EPSG:32607')
detachment13 = gpd.read_file('/Users/mistral/Documents/CUBoulder/Science/Sulzer/VolumeCalculations/MainLoss2013.shp')
detachment13 = detachment13.to_crs('EPSG:32607')

x_limits = [469000,476600]
y_limits = [6833500,6843750]

f,ax = plt.subplots(figsize = (5,6))
ax.imshow(hillshade, extent = dem_extent, cmap = 'Greys')
outline2015.plot(ax = ax, facecolor = 'lightgrey', edgecolor = 'k', linewidth = 0.8)
hummock_area.plot(ax = ax, facecolor = 'pink', linewidth = 0.8)
levee2013.plot(ax = ax, facecolor = 'mediumpurple', linewidth = 0.8, alpha = 0.8)
river2019.plot(ax = ax, facecolor = 'steelblue', linewidth = 0.8)
river2016.plot(ax = ax, facecolor = 'lightskyblue', linewidth = 0.8)
detachment16.plot(ax = ax, facecolor = 'darkred')
detachment15.plot(ax = ax, facecolor = 'darkred', linewidth = 0.8, alpha = 0.9)
detachment13.plot(ax = ax, facecolor = 'orangered', linewidth = 0.8, alpha = 0.9)
detachment15.exterior.plot(ax = ax, color = 'darkred', linewidth = 0.8)
#flowbands.plot(ax = ax, color = 'blue', linewidth = 0.8)
levees.plot(ax = ax, color = 'red', linewidth = 0.8)
ax.set_xlim(x_limits)
ax.set_ylim(y_limits)
f.savefig('zones_overview.pdf')
#f.show()
