import rasterio
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('../ISCE')
import insarhelpers
import geopandas as gpd
import glob
import os
from shapely.geometry import Polygon


band_order = [3,1,0]
res = 5
crs = 'EPSG:32607'

def make_bbox_poly(bounds):
    bbox = Polygon([(bounds[0],bounds[1]),
            (bounds[2],bounds[1]),
            (bounds[2],bounds[3]),
            (bounds[0],bounds[3]),
            (bounds[0],bounds[1])])
    return bbox

def offset_extent(bounds):
    o_e = [bounds[0]-bounds[0],
        bounds[2]-bounds[0],
        bounds[1]-bounds[1],
        bounds[3]-bounds[1]]
    return o_e

# #Planet
# fp = ('/Users/mistral/Documents/CUBoulder/Science/Sulzer+/data/planet/20150813/files')
# imgs = glob.glob(os.path.join(fp,'*_SR_clip.tif'))
#
# def make_img_list(fps):
#     newlist = []
#     for fp in fps:
#         src = rasterio.open(fp)
#         newlist.append(src)
#     return newlist
#
# mosaic_list = make_img_list(imgs)
# planet_rgb, planet_extent = insarhelpers.make_MS_img(mosaic_list, band_order, res = res, brightness_f = 3, contrast_f = 2)

#planet headwall:
src = rasterio.open('/Users/mistral/Documents/CUBoulder/Science/Sulzer/data/GlacierTimelapse/767513_2016-07-29_RE4_3A_Visual_clip.tif')
planet = src.read()
planet_extent = [src.bounds.left, src.bounds.right, src.bounds.bottom, src.bounds.top]

red = insarhelpers.normalize_rgb(planet[0])
green = insarhelpers.normalize_rgb(planet[1])
blue = insarhelpers.normalize_rgb(planet[2])

planet_rgb = (np.dstack((red, green, blue))* 255.0) .astype(np.uint8)

#DEM pre
dem_src_pre = rasterio.open('/Users/mistral/Documents/CUBoulder/Science/Sulzer+/QGIS_Analysis/2016arcticdem_minus2012.tif')
diff_pre = dem_src_pre.read(1)
diff_pre[diff_pre < -1000] = np.nan
pre_extent = [dem_src_pre.bounds.left, dem_src_pre.bounds.right, dem_src_pre.bounds.bottom, dem_src_pre.bounds.top]

#DEM post
dem_src_post = rasterio.open('/Users/mistral/Documents/CUBoulder/Science/Sulzer+/QGIS_Analysis/2019minus2012.tif')
diff_post = dem_src_post.read(1)
diff_post[diff_post < -1000] = np.nan
post_extent = [dem_src_post.bounds.left, dem_src_post.bounds.right, dem_src_post.bounds.bottom, dem_src_post.bounds.top]

#detachment outline
detachment = gpd.read_file('/Users/mistral/Documents/CUBoulder/Science/Sulzer+/QGIS_Analysis/2016_detachment.geojson')

#subplot bounds
box = (469100, 6834200, 470200, 6835400)
zoom_box = gpd.GeoDataFrame([1], geometry = [make_bbox_poly(box)], crs = crs)

#Plot figure
f, (ax1, ax2) = plt.subplots(1,2,figsize = (9, 4))
vis = ax1.imshow(
            planet_rgb,
            extent = planet_extent
)
zoom_box.exterior.plot(ax = ax1, color = 'k', linewidth = 0.8)
ax1.set_xlim([469080.0, 472450.0])
ax1.set_ylim([6833200.0, 6836500.0])
ax1.set_xlabel('Easting (km)')
ax1.set_ylabel('Northing (km)')

vis = ax2.imshow(
            planet_rgb,
            extent = planet_extent
)
detachment.exterior.plot(ax = ax2, color = 'k', linewidth = 0.8)
ax2.set_xlim([469100, 470200])
ax2.set_ylim([6834200, 6835400])
ax2.set_xlabel('X (m)')
ax2.set_ylabel('Y (m)')
f.tight_layout()
#f.show()
f.savefig('vis_2016.pdf')

#DEM Differences
f, (ax1, ax2) = plt.subplots(1,2,figsize = (9, 4))
dem = ax1.imshow(
            diff_pre,
            vmin = -40, vmax = 40,
            cmap = 'seismic_r',
            extent = pre_extent
)
detachment.exterior.plot(ax = ax1, color = 'k', linewidth = 0.8)
ax1.set_xlim([469100, 470200])
ax1.set_ylim([6834200, 6835400])
ax1.set_xlabel('X (m)')
ax1.set_ylabel('Y (m)')

dem = ax2.imshow(
            diff_post,
            vmin = -40, vmax = 40,
            cmap = 'seismic_r',
            extent = post_extent
)
detachment.exterior.plot(ax = ax2, color = 'k', linewidth = 0.8)
ax2.set_xlim([469100, 470200])
ax2.set_ylim([6834200, 6835400])
ax2.set_xlabel('X (m)')
ax2.set_ylabel('Y (m)')
#f.show()
f.savefig('2016_demdiff.pdf', bbox_inches='tight')
