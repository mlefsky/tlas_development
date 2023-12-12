# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
from dask.utils import SerializableLock
import numpy as np
import matplotlib.pyplot as plt
import laspy
import pickle
import os
import subprocess
#import json
#import pdal
import rasterio
#from rasterio.profiles import DefaultGTiffProfile
#import subprocess
import rioxarray
#from rasterio.crs import CRS
#import fiona

def clip(src_fname,trgt_fname,out_fname):

    src  = rioxarray.open_rasterio(src_fname)
    trgt = rioxarray.open_rasterio(trgt_fname)

    xvals=trgt.x.values
    yvals=trgt.y.values
    xmin=min(xvals)
    ymin=min(yvals)
    xmax=max(xvals)
    ymax=max(yvals)
    geometries = [
        {
            'type': 'Polygon',
            'coordinates': [[
                [xmin, ymin],
                [xmin, ymax],
                [xmax, ymax],
                [xmax, ymin],
                [xmin, ymin]
            ]]
        }]

    clipped = src.rio.clip(geometries)
    return(clipped)

    
def mk_pdtm(filename,dtm_filename,out_filename,prefix,**kwargs):
    filename=prefix+filename
    dtm_filename=prefix+dtm_filename
    out_filename=prefix+out_filename
    
    print(filename)
    print(dtm_filename)
    print(out_filename)
    print("prefix: ",prefix)

#                "spatialreference": "EPSG:2193"
#                "returns":"last"

    cmds="""
{
            "pipeline": [
            {
                "type": "readers.las",
                "filename": "$1"
            },
            {
                "type":"filters.assign",
                "assignment":"Classification[:]=0"
            },
            {
                "type":"filters.elm"
            },
            {
                "type":"filters.outlier",
                "method":"statistical",
                "mean_k":8,
                "multiplier":3.0
            },
            {
                "type":"filters.smrf",
                "window":18,
                "threshold":0.30,
                "scalar":0.95
            },
            {
                "type":"filters.range",
                "limits":"Classification[2:2]"
            },
            {
                "type": "writers.gdal",
                "filename": "$3",
                "output_type": "mean",
                "gdaldriver": "GTiff",
                "resolution": 1, 
                "nodata":"-9999",
                "radius": 1,
                "data_type": "float32"
            }

            ]
    }"""
#                "slope":.15,

    ofilename=filename.replace(".las","_pdtm.tif").replace(".laz","_pdtm.tif")
    print("ofilename",ofilename)
    cmds=cmds.replace("$1",filename)
    cmds=cmds.replace("$2",dtm_filename)
    cmds=cmds.replace("$3",ofilename)
    print(cmds)
#    cmds=cmds.replace("$4",out_filename)
    cmds=cmds.replace("\\","\\\\")
    print(">>>>",ofilename)
#    print(">>>>",out_filename)
    print(">>>>>",prefix+'mk_dtm_values.json')
    with open(prefix+'mk_dtm_values.json', 'w') as f:
      f.write(cmds)

    print()
    print("pdal pipeline "+prefix+"mk_dtm_values.json")
    cmd=["pdal", "pipeline",prefix+"mk_dtm_values.json"] 
    print(" ".join(cmd))
    print("past Popen")
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#    subprocess.run("pdal pipeline "+prefix+"mk_dtm_values.json",shell=True)
    print("done")

    return ofilename
    
def mk_pdsm(filename,dsm_filename,out_filename,prefix,**kwargs):
    filename=prefix+filename
    dsm_filename=prefix+dsm_filename
    out_filename=prefix+out_filename
    
    cmds="""
{

    "pipeline":[
        {
            "type": "readers.las",
            "filename": "$1"},
        {
            "type":"filters.range",
            "limits":"returnnumber[0:1]"
        },

        {
            "type": "writers.gdal",
            "filename":"$2",
            "output_type":"max",
            "gdaldriver":"GTiff",
            "resolution": 1,
            "radius": 1,
            "data_type": "float32",
            "nodata":"-9999"

        }
    ]
}
"""
#            "extra_dims":"all",

    ofilename=filename.replace(".las","_pdsm.tif").replace(".laz","_pdsm.tif")
    print(ofilename)
    cmds=cmds.replace("$1",filename)
    cmds=cmds.replace("$2",dsm_filename)
    cmds=cmds.replace("$3",ofilename)
    cmds=cmds.replace("$4",out_filename)
    cmds=cmds.replace("\\","\\\\")

    with open(prefix+'mk_dsm_values.json', 'w') as f:
      f.write(cmds)

    print()
    print("pdal pipeline "+prefix+"mk_dsm_values.json")
    
    # Create a subprocess that has its stdout and stderr streams connected to pipes

    if (1==0):
        print("at Popen")
        cmd=["pdal", "pipeline",prefix+"mk_dsm_values.json"]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print(" ".join(cmd))
        print("past Popen")
        
        # Read the output of the subprocess and get the error message
        output, error = proc.communicate()    
        # Check the output for an error message
        if error:
          print(error)
          
    print("done")

    return ofilename

def combine_dtmdsm(dtm_fname,dsm_fname,prefix):
    
#    dtm  = rioxarray.open_rasterio(prefix+dtm_fname)
#    print(dtm)
#    dsm  = rioxarray.open_rasterio(prefix+dsm_fname)
#    dsm= dtm.squeeze("max")
#    dsm= dsm[1]
#    dtmm= dtm[1]
#    chm_fname=dsm_fname.replace("dsm","chm")
#    print(chm_fname)
#    chm.rio.to_raster(prefix+chm_fname, driver="GTIFF", tiled=True, lock=SerializableLock())  
#    chm = rioxarray.open_rasterio(prefix+chm_fname)
    ofilename=dtm_fname.replace('pdtm','pchm')
    cmd=['gdal_calc.py','-A '+dtm_fname,'-B '+dsm_fname,' --overwrite --calc="B-A" ',' --outfile '+'chm.tif']
    print(" ".join(cmd))  
    subprocess.Popen(cmd)


#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================
import time
#start=time.time()

prefix="/home/lefsky/time_trials/"
dtm_filename="tile_67_136_dtm.tif"
filename="tile_67_136.laz"
dtm_out_filename="tile_67_136_pdtm.tif"
dsm_filename="tile_67_136_dsm.tif"
dsm_out_filename="tile_67_136_pdsm.tif" 

dtm_filename_out=mk_pdtm(filename,dtm_filename,dtm_out_filename,prefix,do_print=True)
dsm_filename_out=mk_pdsm(filename,dsm_filename,dsm_out_filename,prefix,do_print=True)

#print(time.time()-start)
#dsm_filename=mk_pdsm(filename,dsm_filename,dsm_out_filename,prefix,do_print=True)
#tmp=combine_dtmdsm("tile_67_136_pdtm.tif","tile_67_136_pdsm.tif",prefix)

#print(dtm_filename)
#print(dsm_filename)

#==============================================================================
#==============================================================================
#==============================================================================

#=============================================================================
#============================================================================
#=============================================================================

def lasadd_hadtm(filename,dtm_filename,out_filename,prefix,**kwargs):
    filename=prefix+filename
    dtm_filename=prefix+dtm_filename
    out_filename=prefix+out_filename
    print(filename)
    print(dtm_filename)
    print(out_filename)
    print("prefix: ",prefix)
    cmds="""
[
reported.
    {
        "type":"readers.las",
        "filename":"$1",
        "use_eb_vlr":True,
        "compresssion":"laszip"
    },
    {
        "type":"filters.hag_dem",
        "raster":"$2"
    },
    {
        "type": "filters.assign",
        "value": [
            "Z = HeightAboveGround WHERE HeightAboveGround > 0",
            "Z = 0 WHERE Z > 328",
            "Z = 0 WHERE X <= 0"
        ]
    },
    {
        "type":"writers.las",
        "extra_dims":"all",
        "filename":"$3"
    }]

"""

#    {
#        "filename":"$4",
#        "gdaldriver":"la",
#        "output_type":"all",
#        "resolution":"1.0",
#        "type": "writers.gdal"
#    }
#///

    ofilename=prefix+filename.replace(".las","_dtm.las").replace(".laz","_hadtm.laz")
    cmds=cmds.replace("$1",prefix+filename)
    cmds=cmds.replace("$2",prefix+dtm_filename)
    cmds=cmds.replace("$3",ofilename)
    cmds=cmds.replace("$4",prefix+out_filename)
    cmds=cmds.replace("\\","\\\\")

    print(cmds)
    print(">>>>>",prefix+'extract_dtm_values.json')
    with open(prefix+'extract_dtm_values.json', 'w') as f:
      f.write(cmds)

    print()
    print("pdal pipeline "+prefix+"extract_dtm_values.json")
    subprocess.run("pdal pipeline "+prefix+"extract_dtm_values.json",shell=True)
    print("done")
    
    return ofilename



    
def lasvoxelize(filename):
    out_filename=filename.replace(".las","_voxels.las").replace(".laz","_voxels.laz")
    cmd="c:\\LAStools\\bin\\lasvoxel64.exe "+filename+" -empty_voxels -step_xy 14 -step_z 1 -o "+out_filename
    print("cmd: ",cmd)
    subprocess.run(cmd)
    print("done")
    return out_filename

def import_voxels(filename,**kwargs):
    print(kwargs)
    with laspy.open(filename) as fh:
        las = fh.read()
        out_filename=filename.replace(".las","_npv.pkl").replace(".laz","_npv.pkl")
        print(out_filename)

        print('Points from data:', len(las.points))
        print("Stats Calculation")

        x=las.x
        y=las.y
        z=las.z

#        z=las.true_z
        i=las.intensity

        minx=min(x)
        miny=min(y)
        minz=min(z)
        maxx=max(x)
        maxy=max(y)
        maxz=max(z)

        print("x: ",minx,maxx)
        print("y: ",miny,maxy)
        print("z: ",minz,maxz)

        x_index=((x-minx)/14).astype("int")
        y_index=((y-miny)/14).astype("int")
        z_index=(z-minz).astype("int")

        minx=min(x_index)
        miny=min(y_index)
        minz=min(z_index)
        maxx=max(x_index)
        maxy=max(y_index)
        maxz=max(z_index)

        print(minx,maxx)
        print(miny,maxy)
        print(minz,maxz)

        voxels=np.zeros((maxx+1,maxy+1,maxz+1))
        fvoxels=voxels
        print("Filling voxels")

        voxels[x_index.astype(int), y_index.astype(int), z_index.astype(int)] = i.astype(float)

        counts=np.sum(voxels,axis=2)
        voxels[np.isnan(voxels)]=0
        voxels[np.isinf(voxels)]=0
        ocounts=counts
        heights=counts
        print(np.mean(counts))
        print(np.mean(voxels))
        
        for z in range(maxz):
            tmp=voxels[:,:,z]/counts
            tmp[np.isnan(tmp)]=0
            tmp[np.isinf(tmp)]=0
            fvoxels[:,:,z]=tmp
            print(counts.shape)
            print(np.mean(voxels[:,:,z]),np.mean(counts),np.mean(tmp[tmp>0]))
            heights[fvoxels[:,:,z] > 0]=z
            print("sumrat: ",np.sum(counts)/np.sum(ocounts))
        
        fvoxels[np.isnan(fvoxels)]=0
        fvoxels[np.isinf(fvoxels)]=0
        print(np.mean(fvoxels))

        out={'xrange':[minx,maxx],'yrange':[miny,maxy],'range':[minz,maxz],
             "voxels":voxels,\
             "fvoxels":fvoxels,"fname":filename,"ofname":out_filename,\
             "counts":counts,"hmax":heights}
        
        if "save" in kwargs:
            print("Saving...")
            with open(out_filename, 'wb') as f:
                pickle.dump(out, f)
 
        return out

def export_geotiff(image,filename,upperleft,res,crs):
    # Register GDAL format drivers and configuration options with a
    # context manager.
    with rasterio.Env():

        # Write an array as a raster band to a new 8-bit file. For
        # the new file's profile, we start with the profile of the source
        profile = rasterio.profiles.DefaultGTiffProfile({\
        'driver': 'GTiff',
        'interleave': 'band',
        'tiled': True,
        'blockxsize': 256,
        'blockysize': 256,
        'compress': 'lzw',
        'nodata': 0,
        'dtype': type(image[0]),
        'transform' : rasterio.Affine(res, 0.0, upperleft[0],0.0, -1*res, upperleft[1])})
        
        # And then change the band count to 1, set the
        # dtype to uint8, and specify LZW compression.

        profile.update(
            dtype=rasterio.float,
            count=1,
            compress='lzw')
    
#        with rasterio.open(filename, 'w', crs='EPSG:3857', **profile) as dst:
#            dst.write(image.astype(rasterio.float), 1)

def restore_pkl(fname):
    with open(fname, 'rb') as f:
      data = pickle.load(f)
    return data

#filename='F:\\jeffco_collect\\lefsky_usfs_demo_work\\bailey.laz'
#dtm_filename="F:\\jeffco_collect\\lefsky_usfs_demo_work\\bailey_25_dtm.tif"
print("start")
#filename="F:\\jeffco_collect\\laz\\tile_67_136.laz"
#dtm_filename="F:\\jeffco_collect\\dtm\\tile_67_137.tif"
#prefix="/home/lefsky/time_trials/try2/"
#filename="tile_67_136.laz"
#dtm_filename="tile_67_13_dtm.tif"


# -----------------------------------
# ADD Height above ground to laz file
# -----------------------------------

#import time
#start=time.time()
#tmp= add_hag(filename,dtm_filename,prefix,do_print=True)
# 400


# Creates bailey_hag.laz

# -----------------------------------
# CREATE Intermediate voxelized LAZ
# -----------------------------------

#ofilename=voxelize_las("F:\\jeffco_collect\\lefsky_usfs_demo_work\\tile_66_136_z.laz")
#print(ofilename)
#bailey_data=lasvoxelize('F:\\jeffco_collect\\laz\\tile_66_136_hag.laz')#,topo_scale=1,save=True)
#97.9626054763794





# --------------------------------
# IMPORT Intermediate voxelize LAZ
# and
# CALCULATE VOXEL INDICES
# --------------------------------
# If save=TRUE then creates bailey_hag_voxels_npv.pkl
# Returns dict of data products

#baily_np=import_voxels('F:\\jeffco_collect\\laz\\tile_67_136_hag_voxels.laz',save=True)
#print("finish:",time.time()-start)

# --------------------------------
# RESTORE pkl file and show elements of data dict
#data=restore_pkl("bailey_hag_voxels_npv.pkl")
#data.keys()
#dict_keys(['xrange', 'yrange', 'range', 'voxels', 'fvoxels', 'fname', 'ofname', 'counts', 'hmax'])







# OBSOLETE
# OBSOLETE
# OBSOLETE
# voxels=voxelize("bailey_hag.laz",save=True)
##
##def voxelize(filename,**kwargs):
##    print(kwargs)
##    with laspy.open(filename) as fh:
##        las = fh.read()
##        out_filename=filename.replace(".las","_npv.pkl").replace(".laz","_npv.pkl")
##        print(out_filename)
##
##        print(las)
##        print('Points from data:', len(las.points))
##        print("Stats Calculation")
##
##        x=las.x
##        y=las.y
##        z=las.z
##        #z=las.true_z
##        #i=las.intensity
##
##        minx=min(x)
##        miny=min(y)
##        minz=min(z)
##        maxx=max(x)
##        maxy=max(y)
##        maxz=max(z)
##
##        print("x: ",minx,maxx)
##        print("y: ",miny,maxy)
##        print("z: ",minz,maxz)
##
##        if ("res" not in kwargs):
##            res=14
##            
##        x_index=((x-minx)/res).astype("int")
##        y_index=((y-miny)/res).astype("int")
##        z_index=(z-minz).astype("int")
##
##
##        minx=min(x_index)
##        miny=min(y_index)
##        minz=min(z_index)
##        maxx=max(x_index)
##        maxy=max(y_index)
##        maxz=max(z_index)
##        
##        voxels=np.histogramdd(np.array([x_index, y_index, z_index]), range=[(0,maxx),(0,maxy),(0,maxz)])#, weights=None)
##        print(voxels.shape)
##        fvoxels=voxels
##
##        print(minx,maxx)
##        print(miny,maxy)
##        print(minz,maxz)
##
###        voxels=np.zeros((maxx+1,maxy+1,maxz+1))
###        voxels1=voxels
###        fvoxels=voxels
###        print("Filling voxels")
##
###        voxels[x_index.astype(int), y_index.astype(int), z_index.astype(int)] = i.astype(float)
##    
##        counts=np.sum(voxels,axis=2).shape
##        
##        for z in range(maxz):
##            fvoxels[:,:,z]=voxels[:,:,z]/counts
##
##        out={'xrange':[minx,maxx],'yrange':[miny,maxy],'range':[minz,maxz],"voxels":voxels,\
##             "fvoxels":fvoxels,"fname":filename,"ofname":out_filename,\
##             "counts":counts}
##        
##        if "save" in kwargs:
##            print("Saving...")
##            with open(out_filename, 'wb') as f:
##                pickle.dump(out, f)
##    
##        return voxels
##
##
