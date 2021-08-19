import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.geometry import Polygon
import earthpy
import earthpy.plot
import earthpy.spatial
import rasterio
import rasterio.merge
import rasterio.plot
import numpy as np
import glob
import os
from PIL import Image, ImageEnhance

# import geol map
bbox = (630574, 1345138, 671174, 1371206)
geology = gpd.read_file('/Users/mistral/Documents/CUBoulder/Science/spatial_base_data/AKgeol_web_shp/AKStategeol_poly.shp', bbox = bbox)
geology = geology.to_crs('EPSG:32607')

#import DEM
dem = '/Users/mistral/Documents/CUBoulder/Science/Sulzer+/data/DEM/regionalIFSAR_UTM_10m.tif'
with rasterio.open(dem) as src:
    elevation = src.read(1)
    elevation[elevation < 0] = np.nan

dem_extent = [373271.726, 521821.726, 6743453.939, 6894553.939]
hillshade = earthpy.spatial.hillshade(elevation)

# Optical image (post 2015 flow)
im_path = '/Users/mistral/Documents/CUBoulder/Science/Sulzer+/data/planet/20160731/files'
im_fp = glob.glob(os.path.join(im_path,'*Analytic_SR_clip.tif'))
mosaic_list = []

#way to scale all planet images in same way for better plotting...
src1 = rasterio.open(im_fp[0])

for fp in im_fp:
    src = rasterio.open(fp)
    mosaic_list.append(src)

mosaic, out_transform = rasterio.merge.merge(mosaic_list, res = (5,5))
mosaic_extent = [out_transform[2], out_transform[2]+(mosaic.shape[2]*out_transform[0]),
                out_transform[5]+(mosaic.shape[1]*out_transform[4]), out_transform[5]]

def normalize(array):
    """Normalizes numpy arrays into scale 0.0 - 1.0"""
    array_min, array_max = array.min(), array.max()
    return ((array - array_min)/(array_max - array_min))

red = normalize(mosaic[0])
green = normalize(mosaic[1])
blue = normalize(mosaic[2])
nir = normalize(mosaic[3])

planet_rgb = (np.dstack((nir, green, red))* 255.0) .astype(np.uint8)
img = Image.fromarray(planet_rgb) # turn in to PIL image
brightness = ImageEnhance.Brightness(img)
b_factor = 2
planet_bright = brightness.enhance(b_factor)
contrast = ImageEnhance.Contrast(planet_bright)
c_factor = 1
planet_contrast = contrast.enhance(c_factor)


# Field site frames:
#Flat Creek
fc_poly = Polygon([(468660,6833110),(472230,6843520),(477920,6843520),(472550,6833110),(468660,6833110)])
fc_poly_gdf = gpd.GeoDataFrame([1], geometry = [fc_poly], crs = 'EPSG:32607')
#West site
ws_poly = Polygon([(475110,6833850),(475110,6835270),(477220,6835270),(477220,6833850),(475110,6833850)])
ws_poly_gdf = gpd.GeoDataFrame([1], geometry = [ws_poly], crs = 'EPSG:32607')

x_limits = [466000,480000]
y_limits = [6831000,6846000]

f,ax = plt.subplots(figsize = (8,6))
ax.imshow(planet_contrast, extent = mosaic_extent, vmin = 0.45, vmax = 0.55)
fc_poly_gdf.boundary.plot(ax=ax, color="k", linewidth = 1)
ws_poly_gdf.boundary.plot(ax=ax, color="k", linewidth = 1)
ax.text(474030,6843650,'Flat Creek', color = 'k')
ax.text(475410,6835400,'East Site', color = 'w')
ax.set_xlim(x_limits)
ax.set_ylim(y_limits)
#f.show()
f.savefig('site_overview_2016_2.pdf', dpi=300)

# Geologic Map on hillshade
f,ax = plt.subplots(figsize = (12,6))
ax.imshow(hillshade, extent = dem_extent, cmap = 'Greys')
geology.plot(ax = ax, column = 'STATE_UNIT', categorical=True, legend=True,
    edgecolor='0.2', linewidth = 0.2,
    legend_kwds={'title': "Geology", 'ncol':1, 'bbox_to_anchor':(2,1.05), 'frameon':False},
    cmap = 'Paired', alpha=0.6)
fc_poly_gdf.boundary.plot(ax=ax, color="k", linewidth = 1)
ws_poly_gdf.boundary.plot(ax=ax, color="k", linewidth = 1)
ax.set_xlim(x_limits)
ax.set_ylim(y_limits)
ax.text(474030,6843650,'Flat Creek', color = 'k')
ax.text(475410,6835400,'East Site', color = 'k')
f.tight_layout()
f.show()
#f.savefig('geology.pdf')

#Geographical overview
#import rgi
rgi = gpd.read_file('/Users/mistral/Documents/CUBoulder/Science/spatial_base_data/01_rgi60_Alaska/01_rgi60_Alaska.shp')
rgi = rgi.to_crs('EPSG:32607')

plot_window = Polygon([(x_limits[0],y_limits[0]),(x_limits[0],y_limits[1]),(x_limits[1],y_limits[1]),(x_limits[1],y_limits[0]),(x_limits[0],y_limits[0])])
plot_window_gdf = gpd.GeoDataFrame([1], geometry = [plot_window], crs = 'EPSG:32607')

f, ax = plt.subplots(figsize = (8,6))
ax.imshow(hillshade, cmap = 'Greys', extent = dem_extent)
rgi.plot(ax = ax, facecolor = 'lightblue', alpha = 0.5)
plot_window_gdf.boundary.plot(ax = ax, color = 'k')
ax.set_xlim([393000,485000])
ax.set_ylim([6790000,6865000])
f.show()
#f.savefig('regional_overview.pdf')
