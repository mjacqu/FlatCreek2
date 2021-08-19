import numpy as py
import sys
sys.path.append('../ISCE')
import insarhelpers
import geopandas as gpd
import matplotlib.pyplot as plt
import glob
import os
from PIL import Image, ImageEnhance
import rasterio
import rasterio.merge
import rasterio.plot
import numpy as np
from shapely.geometry import Polygon
import geopandas as gpd


datapath = '/Users/mistral/Documents/CUBoulder/Science/Sulzer+/data/ChrisLarsen'
imgs16 = glob.glob(os.path.join(datapath,'WhiteRiver_June1_2016_ortho*.tif'))
imgs19 = glob.glob(os.path.join(datapath,'Aug30_2019_Sulzer*.tif')) #flat creek and east site
eastsite_imgs = glob.glob(os.path.join(datapath,'Aug30_2019_Sulzer_east_ortho*.tif'))

def make_img_list(fps):
    newlist = []
    for fp in fps:
        src = rasterio.open(fp)
        newlist.append(src)
    return newlist

mosaic_list16 = make_img_list(imgs16)
mosaic_list19 = make_img_list(imgs19)
eastsite = make_img_list(eastsite_imgs)

def make_bbox_poly(bounds):
    bbox = Polygon([(bounds[0],bounds[1]),
            (bounds[2],bounds[1]),
            (bounds[2],bounds[3]),
            (bounds[0],bounds[3]),
            (bounds[0],bounds[1])])
    return bbox

band_order = [0,1,2]
res = 5
crs = 'EPSG:32607'
bounds16 = (469000,6833500,479500,6844000)
hummocks_bounds = (474791, 6840032, 475066, 6840214)
tree_bounds= (472420, 6842385, 472630, 6842555)
zoom_bounds = (473400,6839400,474400, 6840200)
deflation_bounds = (473944,6839191, 474202, 6839371)
flowlines_bounds = (474074, 6838819, 474297, 6839060)
stripes_bounds = (471765,6836367,472040,6836549)
water_bounds = (473905,6839729,474204,6839968)
scratches_bounds = (471588,6835718,471695,6835812)

hummocks_box = gpd.GeoDataFrame([1], geometry = [make_bbox_poly(hummocks_bounds)], crs = crs)
tree_box = gpd.GeoDataFrame([1], geometry = [make_bbox_poly(tree_bounds)], crs = crs)
zoom_box = gpd.GeoDataFrame([1], geometry = [make_bbox_poly(zoom_bounds)], crs = crs)
deflation_box = gpd.GeoDataFrame([1], geometry = [make_bbox_poly(deflation_bounds)], crs = crs)
flowlines_box = gpd.GeoDataFrame([1], geometry = [make_bbox_poly(flowlines_bounds)], crs = crs)
stripes_box = gpd.GeoDataFrame([1], geometry = [make_bbox_poly(stripes_bounds)], crs = crs)
water_box = gpd.GeoDataFrame([1], geometry = [make_bbox_poly(water_bounds)], crs = crs)


ortho2016, extent = insarhelpers.make_MS_img(mosaic_list16, band_order, res = res, bounds = bounds16)
ortho2019, extent = insarhelpers.make_MS_img(mosaic_list19, band_order, res = res)
ortho_east, east_extent = insarhelpers.make_MS_img(eastsite, band_order, res = res)
hummocks, hummocks_extent = insarhelpers.make_MS_img(mosaic_list19, band_order, bounds = hummocks_bounds)
trees, trees_extent = insarhelpers.make_MS_img(mosaic_list19, band_order, bounds = tree_bounds)
deflation16, deflation_extent = insarhelpers.make_MS_img(mosaic_list16, band_order, bounds = deflation_bounds)
deflation19, deflation_extent = insarhelpers.make_MS_img(mosaic_list19, band_order, bounds = deflation_bounds)
flowlines, fl_extent = insarhelpers.make_MS_img(mosaic_list16, band_order, bounds = flowlines_bounds)
stripes, stripes_extent = insarhelpers.make_MS_img(mosaic_list19, band_order, bounds = stripes_bounds)
waterbodies, water_extent = insarhelpers.make_MS_img(mosaic_list16, band_order, bounds = water_bounds)
scratches, scratches_extent = insarhelpers.make_MS_img(mosaic_list19, band_order, bounds = scratches_bounds)

#load polygons
ertlines = gpd.read_file('/Users/mistral/Documents/CUBoulder/Science/Sulzer+/QGIS_Analysis/ERT_Lines.geojson')
ertlines = ertlines.to_crs('EPSG:32607')

sedsamples = gpd.read_file('/Users/mistral/Documents/CUBoulder/Science/Sulzer+/QGIS_Analysis/SedSamplesFlatCreek.geojson')
sedsamples = sedsamples.to_crs('EPSG:32607')

glaciers = gpd.read_file('/Users/mistral/Documents/CUBoulder/Science/spatial_base_data/01_rgi60_Alaska/01_rgi60_Alaska.shp')
glaciers = glaciers.to_crs('EPSG:32607')


levees = gpd.read_file('/Users/mistral/Documents/CUBoulder/Science/Sulzer+/QGIS_Analysis/levees.geojson')
water = gpd.read_file('/Users/mistral/Documents/CUBoulder/Science/Sulzer+/QGIS_Analysis/water.geojson')
flowbands = gpd.read_file('/Users/mistral/Documents/CUBoulder/Science/Sulzer+/QGIS_Analysis/pressure_ridges.geojson')
runout13 = gpd.read_file('/Users/mistral/Documents/CUBoulder/Science/Sulzer+/QGIS_Analysis/2013_runout.geojson')
runout15 = gpd.read_file('/Users/mistral/Documents/CUBoulder/Science/Sulzer+/QGIS_Analysis/2015_runout.geojson')
molards = gpd.read_file('/Users/mistral/Documents/CUBoulder/Science/Sulzer+/QGIS_Analysis/molards.geojson')
eastdposit = gpd.read_file('/Users/mistral/Documents/CUBoulder/Science/Sulzer+/QGIS_Analysis/eastdeposit_outline.geojson')


#create mask for ortho image
mask = ((np.array(ortho2019)[..., 0] != 0)*255).astype('uint8')
ortho2019.putalpha(Image.fromarray(mask, mode="L"))

#overview plot
f, ax = plt.subplots(figsize = (7.5, 9))
ax.imshow(ortho2019, extent = extent)
hummocks_box.exterior.plot(ax = ax, color = 'k', linewidth = 0.8)
tree_box.exterior.plot(ax = ax, color = 'k', linewidth = 0.8)
deflation_box.exterior.plot(ax = ax, color = 'k', linewidth = 0.8)
flowlines_box.exterior.plot(ax = ax, color = 'k', linewidth = 0.8)
stripes_box.exterior.plot(ax = ax, color = 'k', linewidth = 0.8)
water_box.exterior.plot(ax = ax, color = 'k', linewidth = 0.8)
runout13.exterior.plot(ax = ax, color = 'yellow', linewidth = 0.8)
runout15.exterior.plot(ax = ax, color = 'red', linewidth = 0.8)
ertlines.plot(ax = ax, color = 'w', linewidth = 1)
sedsamples.plot(ax = ax, color = 'k', markersize = 8)
eastdposit.plot(ax = ax, color = 'w', linewidth = 0.8)
#levees.plot(ax = ax, color = 'gold', linewidth = 0.5)
flowbands.plot(ax = ax, color = 'red', linewidth = 0.5)
water.plot(ax = ax, facecolor = 'deepskyblue', edgecolor = 'deepskyblue', linewidth = 0.5)
#molards.plot(ax = ax, color = 'lime', markersize = 1)
runout13.exterior.plot(ax = ax, color = 'w', linewidth = 0.3)
ax.set_xlim([469000,477000])
ax.set_ylim([6833500,6844000])
ax.ticklabel_format(useOffset=None, style = 'plain')
ax.set_xlabel('Easting (m)', fontsize = 16)
ax.set_ylabel('Northing (m)', fontsize = 16)
f.tight_layout()
#f.show()
f.savefig('ortho3.pdf')

#get m-extents:
def offset_extent(bounds):
    o_e = [bounds[0]-bounds[0],
        bounds[2]-bounds[0],
        bounds[1]-bounds[1],
        bounds[3]-bounds[1]]
    return o_e

#hummocks plot
f,ax = plt.subplots()
ax.imshow(hummocks, extent = offset_extent(hummocks_bounds))
ax.set_xlabel('X (m)')
ax.set_ylabel('Y (m)')
f.tight_layout()
#f.show()
f.savefig('hummocks.pdf')

#trees plot
f,ax = plt.subplots()
ax.imshow(trees, extent = offset_extent(tree_bounds))
ax.set_xlabel('X (m)')
ax.set_ylabel('Y (m)')
f.tight_layout()
f.show()
#f.savefig('trees.pdf')

#deflation plot
f, (ax1, ax2) = plt.subplots(1,2, sharey = True, figsize = (10,5))
ax1.imshow(deflation16, extent = offset_extent(deflation_bounds))
ax2.imshow(deflation19, extent = offset_extent(deflation_bounds))
ax1.set_xlabel('X (m)')
ax2.set_xlabel('X (m)')
ax1.set_ylabel('Y (m)')
f.tight_layout()
f.show()
#f.savefig('deflation.pdf')

#flowbands plot
f,ax = plt.subplots()
ax.imshow(flowlines, extent = offset_extent(flowlines_bounds))
ax.set_xlabel('X (m)')
ax.set_ylabel('Y (m)')
f.tight_layout()
f.show()
#f.savefig('flowbands.pdf')

#stripes plot
f,ax = plt.subplots()
ax.imshow(stripes, extent = offset_extent(stripes_bounds))
ax.set_xlabel('X (m)')
ax.set_ylabel('Y (m)')
f.tight_layout()
#f.show()
f.savefig('stripes.pdf')

#waterbodies plot
f,ax = plt.subplots()
ax.imshow(waterbodies, extent = offset_extent(water_bounds))
ax.set_xlabel('X (m)')
ax.set_ylabel('Y (m)')
f.tight_layout()
f.show()
#f.savefig('waterbodies.pdf')

#scratches plot
f,ax = plt.subplots()
ax.imshow(scratches, extent = offset_extent(scratches_bounds))
ax.set_xlabel('X (m)')
ax.set_ylabel('Y (m)')
f.tight_layout()
#f.show()
f.savefig('scratches.pdf')

# DEM differences

#load DEM difference
dem_src = rasterio.open('/Users/mistral/Documents/CUBoulder/Science/Sulzer+/QGIS_Analysis/2019minus2016_clipped.tif')
dem_diff = dem_src.read(1)
dem_extent = [dem_src.bounds.left, dem_src.bounds.right, dem_src.bounds.bottom, dem_src.bounds.top]

dem_diff = dem_diff - 0.3 #correct by median value over runout zone
dem_diff[dem_diff < -1000] = np.nan

#boxes and close up views
dem_defl = (473900,6839100, 474200, 6839375)
dem_defl_box = gpd.GeoDataFrame([1], geometry = [make_bbbox_poly(dem_defl)], crs = crs)

ert_line = (471650, 6835950, 472000, 6836250)
ert_line_box = gpd.GeoDataFrame([1], geometry = [make_bbbox_poly(ert_line)], crs = crs)

glacier_view = (469700, 6833900, 471875, 6835925)
glacier_view_box = gpd.GeoDataFrame([1], geometry = [make_bbbox_poly(glacier_view)], crs = crs)

f.clf()
plt.close()

f, ax = plt.subplots(figsize = (7.5, 9)) #replot for publication with correct boxes
ax.imshow(ortho2019, extent = extent)
dem = ax.imshow(
            dem_diff,
            vmin = -10, vmax = 10,
            cmap = 'seismic_r',
            extent = dem_extent
)
#ertlines.plot(ax = ax, color = 'k', linewidth = 1)
dem_defl_box.exterior.plot(ax = ax, color = 'k', linewidth = 0.8)
ert_line_box.exterior.plot(ax = ax, color = 'k', linewidth = 0.8)
glacier_view_box.exterior.plot(ax = ax, color = 'k', linewidth = 0.8)
ax.set_xlim([dem_src.bounds.left,dem_src.bounds.right])
ax.set_ylim([dem_src.bounds.bottom, dem_src.bounds.top])
ax.ticklabel_format(useOffset=None, style = 'plain')
ax.set_xlabel('Easting (m)', fontsize = 14)
ax.set_ylabel('Northing (m)', fontsize = 14)
cb = f.colorbar(dem)
cb.ax.minorticks_off()
cb.set_label('Elevation change')
f.tight_layout()
#f.show()
f.savefig('DemDiff.pdf')

#dem deflation plot
f,ax = plt.subplots()
dem = ax.imshow(
        dem_diff,
        vmin = -4, vmax = 4,
        cmap = 'seismic_r',
        extent = dem_extent
)
cb = f.colorbar(dem, orientation = 'horizontal')
cb.ax.minorticks_off()
cb.set_label('Elevation change')
ax.set_xlim([dem_defl[0],dem_defl[2]])
ax.set_ylim([dem_defl[1],dem_defl[3]])
ax.set_xlabel('X (m)')
ax.set_ylabel('Y (m)')
f.tight_layout()
#f.show()
f.savefig('DemDiff_deflation.pdf')

#ert line plot
f,ax = plt.subplots()
ax.imshow(ortho2019, extent = extent)
dem = ax.imshow(
        dem_diff,
        vmin = -8, vmax = 8,
        cmap = 'seismic_r',
        extent = dem_extent
)
ax.set_xlim([ert_line[0],ert_line[2]])
ax.set_ylim([ert_line[1],ert_line[3]])
ax.set_xlabel('X (m)')
ax.set_ylabel('Y (m)')
cb = f.colorbar(dem, orientation = 'horizontal')
cb.ax.minorticks_off()
cb.set_label('Elevation change')
f.tight_layout()
#f.show()
f.savefig('DemDiff_ertline.pdf')

#glacier view plot
f,ax = plt.subplots()
ax.imshow(ortho2019, extent = extent)
dem = ax.imshow(
        dem_diff,
        vmin = -20, vmax = 20,
        cmap = 'seismic_r',
        extent = dem_extent
)
glaciers.exterior.plot(ax = ax, color = 'k')
ax.set_xlim([glacier_view[0],glacier_view[2]])
ax.set_ylim([glacier_view[1],glacier_view[3]])
ax.set_xlabel('X (m)')
ax.set_ylabel('Y (m)')
cb = f.colorbar(dem)
cb.ax.minorticks_off()
cb.set_label('Elevation change')
f.tight_layout()
#f.show()
f.savefig('DemDiff_glacierview.pdf')
