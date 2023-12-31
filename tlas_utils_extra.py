
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
    print(">>>>>",prefix+'extract_dtm.json')
    with open(prefix+'extract_dtm.json', 'w') as f:
      f.write(cmds)

    print()
    print("pdal pipeline "+prefix+"extract_dtm.json")
#    subprocess.run("pdal pipeline "+prefix+"extract_dtm.json",shell=True)
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











