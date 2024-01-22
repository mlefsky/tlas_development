# -*- coding: utf-8 -*-
# """
# Spyder Editor

# This is a temporary script file.
# """
# #!/usr/bin/env python
# import rsgislib
# from rsgislib import imageutils
# # Create list of images
 bandNamesList=['nodup_100_22_density_lt1','nodup_100_22_density_lt1','nodup_100_22_cover',
 'nodup_100_22_pdtm','nodup_100_22_pdsm','nodup_100_22_chm','nodup_512_11_density_lt1',
 'nodup_512_11_density','nodup_512_11_cover','nodup_512_11_pdtm','nodup_512_11_pdsm','nodup_512_11_chm',
 'nodup_669_8_density_lt1','nodup_669_8_density','nodup_669_8_cover','nodup_669_8_pdtm','nodup_669_8_pdsm','nodup_669_8_chm',
 'nodup_985_4_density_lt1','nodup_985_4_density','nodup_985_4_cover','nodup_985_4_pdtm','nodup_985_4_pdsm','nodup_985_4_chm',
 'sub_01_chm','sub_01_cover','sub_01_density_lt1','sub_01_density','sub_01_pdsm','sub_01_pdtm']
print(len(bandNamesList))
# print(len(bandNamesList))

# imageList = ["nodup_100_22_stack.tif",
# "odup_512_11_stack.tif",
# "nodup_669_8_stack.tif",
# "nodup_985_4_stack.tif",
# "tile_66_136_sub_01_stack.tif"]

# # Set output image
# outputImage = 'density_stack.tif'
# # Set format and type
# gdalFormat = 'GTiff'
# dataType = rsgislib.TYPE_32FLOAT
# # Stack
# imageutils.stack_img_bands(imageList, bandNamesList, outputImage, -9999, -9999, gdalFormat, dataType)



from osgeo import gdal
import geopandas
import rasterio
import matplotlib.pyplot as plt
import  numpy
from shapely.geometry import Point

file="/home/ubuntu/density_comparison_v2/multistack2.tif"
src = gdal.Open(file)
ulx, xres, xskew, uly, yskew, yres  = src.GetGeoTransform()
lrx = ulx + (src.RasterXSize * xres)
lry = uly + (src.RasterYSize * yres)
for i in range(1, src.RasterCount + 1):
    band = src.GetRasterBand(1).GetDescription()
    print(i,band)
print(ulx, xres, xskew, uly, yskew, yres,lrx,lry)
print(lrx-ulx,uly-lry)
print(bands)

xout=[]
yout=[]
skip=20

for x in range(int(ulx),int(lrx),skip):
    for y in range(int(lry),int(uly),skip):
        xout.append(x)
        yout.append(y) 
#    print(lrx,x,uly)

# Create sampling points
points = []

for ix,x in enumerate(xout):
    
    points.append(Point(xout[ix],yout[ix]))

lix=list(range(len(xout)))
#print(lix)
gdf = geopandas.GeoDataFrame(lix, geometry=points)
                             
gdf.head()

src = rasterio.open(file)

coord_list = [(x, y) for x, y in zip(gdf["geometry"].x, gdf["geometry"].y)]
gdf["value"] = [x for x in src.sample(coord_list)]
print(x)
gdf.head()
numpy.set_printoptions(linewidth=100000)
gdf.to_csv("/home/ubuntu/density_comparison_v2/multistack_out.csv")
