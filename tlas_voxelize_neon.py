import sys
sys.path.append("/home/lefsky/time_trials/Lidar-Notebooks-main/bin")
#    sys.path.append('../../bin/')
from LabLidar_Functions import calccover, calcPercentileHeights, canopyLayerMetrics
from LabLidar_Classes import Cloud
import geopandas as gpd
import pandas as pd
import numpy as np
import concurrent.futures
from pathlib import Path
import laspy
import time
import matplotlib.pyplot as plt
import seaborn as sns
import pickle
from shapely.geometry import Polygon
from scipy.ndimage import gaussian_filter1d
from scipy.stats import gaussian_kde
from scipy.signal import find_peaks, peak_widths
from scipy.interpolate import interp1d
from glob import glob
import os

def flatten_comprehension(matrix):
    return [item for row in matrix for item in row]   


#"https://github.com/pbb2291/Lidar-Notebooks"




def tlas_voxelize_init(lasinputs,las_directory,output_directory,max_height,xysize,units="m"):
    # Import Dependencies
    
    # makes matplotlib plots big
    plt.rcParams['figure.figsize'] = [8, 6]
    plt.rcParams.update({'font.size': 14})
    
    # # # USER INPUTS
    
    # Input directory of las files to compute metrics over
#    ld = Path('/n/davies_lab/Lab/LabLidarScripts/data/out/test/nkhulu/las_clipped/')
     
   # make las inputs
    if lasinputs=="":
        lasinputs=[]
        lasinputs.append(glob(las_directory+'*.las'))
        lasinputs.append(glob(las_directory+'*.laz'))
        lasinputs=flatten_comprehension(lasinputs)
        print(lasinputs)
            
    od_metrics = Path(output_directory)
    
    # EPSG code of the las files, as a string
    # Kruger is 32736 (WGS84 UTM 36S)
    # Mpala is 32637 (WGS84 UTM 37N)
    # Selenkay is 32737 (WGS84 UTM37S)
#    epsg='32736'
    
    # Max height of voxel stacks 
    # NOTE: Set this to be just above max height of your trees in meters.
    stackheight = 20
    
    
    # Horizontal Res of Grid (XY pixel size)
    xysize = xysize
    
    # Vertical step size for metrics
    # NOTE: This defines the vertical bin size in meters (how "thick" each voxel is).
    verticalres = 1
    
    # Set the ground threshold in meters (i.e. below this height treat points as ground).
    # # # 
    # EXPLANATION:
    # You can use groundthreshold to account for errors in relative accuracy
    # For example, if the rel. accuracy of ground is about 0.06 m (6 cm) between flightlines,
    # You could add in a voxel bin that extends from 0-0.06 m, 
    # treating all points with a height in that range as ground points.
    # (so any hit below 0.06 m counts as ground).
    # If you prefer to use all points above 0 m, just set groundthres to 0.
    # Note: groundthreshold can also work with negative heights,
    # setting it to -0.05 for instance would treat points 
    # with negative height values (between 0 and -5cm) as ground points.
    # # # 
    groundthreshold = 0.25
    
    # height col
    heightcol = 'HeightAboveGround'
    
    # Set Complexity Metric Parameters
    # Note: Please leave these as their default values,
    # unless you have a strong reason to change them. 
    # They have been parameterized for savanna environments. 
    
    # Set method for calculation of peaks/layers (options = 'kde' or 'gauss1d')
    method = 'gauss1d'
    
    # Set smoothing sigma if using 'gauss1d"
    sigma=0.5
    
    # set relative height for top of herb layer calc in canopyLayerMetrics
    rh = 0.9
    
    # # # END USER INPUTS
    #====================================================================
    for li in lasinputs:
        print(lasinputs)
        if not os.path.exists(li):
            print(f'Warning: {li} does not exist')
    
    if len(lasinputs) < 1:
        print(f'Warning: Empty input directory- no las files found')
        
    if not od_metrics.exists():
        print('Warning: Output directory does not exist')
    
    #====================================================================
    # Start Voxelizing
    
    # Make voxel height bins
    # Calc Cover for height bins
    nbins = ((stackheight - 0) / verticalres) + 1
    heightbins = np.linspace(0, stackheight, int(nbins))
    
    if groundthreshold > 0:
        # insert the groundthres into the array (right above 0)
        heightbins = np.insert(heightbins, 1, groundthreshold)
    if groundthreshold < 0:
        # insert the groundthres into the array (right below 0)
        heightbins = np.insert(heightbins, 0, groundthreshold)
    
    
    outstruct={"lasinputs":lasinputs,
        "las_directory":las_directory,
        "od_metrics":od_metrics,
        "stackheight": stackheight,
        "xysize":xysize,
        "verticalres":1,
        "groundthreshold":0.25,
        "heightcol":'HeightAboveGround',
        "method" : "gauss1d",
        "sigma":0.5,
        "rh":0.9,
        "nbins":nbins,
        "heightbins":heightbins}
#        "indices":indices}
        
    return(outstruct)
    #===================================================================
    # Wrapper function for using parallel processing and calccover function 
    # Notice that is calls lc as the first argument
    # need to write it this way in order to use concurrent futures parallel processing below

    


def calccover_parallel(lc,index, struct):

    # make a True/False array 
    # for all points within the current grid cell
#    return(lc,index,struct)
    idx_bool = lc.grid_dict['idx_points'] == lc.grid_dict['idx_cells'][index]
    return(idx_bool,lc,index)
    # Subset Points
    p = lc.las.points[idx_bool]

    # Get height array
    # Note: this is slightly different from the "heights" output below
    h = struct['heightcol']
    
    # Remove high noise points above the canopy
    h = h[h<=struct['stackheight']]
#    las_directory,output_directory,max_height,xysize,h_units="m"):
    # Import Dependencies
    try:

        # Calculate Cover
        cover = calccover(points=p,
                          heightbins=struct['heightbins'], 
                          step=struct['verticalres'],
                          groundthres=struct['groundthreshold'],
                          heightcol=struct['heightcol'],
                          hmax=struct['max_height'])
        
    except Exception as e:

        print(f"Cover Calc. - {e.__class__} for {lc.lasf}: \n")
        print(f"\t{e}\n")

    try: 
        
        # Calculate height statistics, and return an array o"cover":cover, "perc", perc, "heights":heights, "complexity":complexity}f the point heights above groundthreshold
        perc, heights = calcPercentileHeights(points=p,
                                              groundthres=struct['groundthreshold'],
                                              returnHeights=True,
                                              heightcol=struct['heightcol'],
                                              hmax=struct['stackheight'])
    except Exception as e:

        print(f"Percentile Calc. - {e.__class__} for {lc.lasf}: \n")
        print(f"\t{e}\n")

    try:
        
        # Compute complexity metrics
        complexity = canopyLayerMetrics(h=h,
                                        hbins=struct['heightbins'],
                                        method=struct['method'],
                                        smoothsigma=struct['sigma'],
                                        rel_height=struct['rh'],
                                        groundthreshold=struct['groundthreshold'])
    except Exception as e:

        print(f"Complexity Calc. - {e.__class__} for {lc.lasf}: \n")
        print(f"\t{e}\n") 

    # Return cover dict, percentile dict, height array, and complexity metrics
    return {"cover":cover, "perc": perc, "heights":heights, "complexity":complexity}


def tlas_voxelize(lasinputs,las_directory,output_directory,max_height,xysize,units="m"):
       
    params=tlas_voxelize_init(lasinputs,las_directory,output_directory,max_height,xysize,units="m")
    
    
    for lasf in lasinputs:    
        
        ### STEP 1: Load in Cloud 
        startproj = time.time()
    
        # Make a las cloud class, and grid it
        lc = Cloud(lasf=lasf,
                   gridsize=params["xysize"],
                   vsize=params["verticalres"],
                   heightcol=params["heightcol"],
                   maxh=params["stackheight"])
    
        # get project string from file name for saving below
#        projstr = Path(lasf).name.split('.')[0]
    
        end = time.time()

        projstr=""        
        print(f'Loaded {projstr} cloud, time elapsed: {end - startproj}\n')
        
        ### STEP 2: Make the grid
        start = time.time()
        
        lc.makegrid()
    
        end = time.time()
        
        print(f'Grid created, time elapsed: {end - start}\n')
        
        ### STEP 2: Compute Cover, FHP, and Percentiles Metrics Over the Cloud's Grid
        start = time.time()
    
        # initialize dictionaries for output 
        lc.cover_dict = {}
        lc.perc_dict = {}
        lc.height_dict = {}
        lc.complexity_dict = {}
    
        # set the cell indices to loop over in parallel
        indices = lc.grid_dict['idx_cells']
        idx_bool,lc,index=calccover_parallel(lc,indices, params)
        return (idx_bool,lc,index)
        numcells = len(lc.grid_dict['idx_cells']) 
        
        print(f'Starting parallel processing of {numcells} pixels.\n')
        # Assuming 2 seconds for each run, and 40 processes as once (should be a conservative estimate)
        print(f'\tApproximate time to completion of {lasf.name}:{np.round(((numcells/40)*2)/36000, 3)} hrs')
    
        ## Use concurrent futures to compute cover over each cell in parallel
        with concurrent.futures.ProcessPoolExecutor(max_workers=None) as executor:
            for cphc, x, y in zip(executor.map(calccover_parallel, indices),
                                     lc.grid_dict['x_cells'],
                                     lc.grid_dict['y_cells']):
    
                try:
    
                    # Stick the cover, perc, and heights inside the metrics dictionary
                    # with x and y location as tuple keys
                    lc.cover_dict[(x, y)] = cphc[0]
                    lc.perc_dict[(x, y)] = cphc[1]
                    lc.height_dict[(x, y)] = np.round(cphc[2], decimals=3)
                    lc.complexity_dict[(x, y)] = cphc[3]
    
                except Exception as e:
    
                    print(f"Saving metrics error - {e.__class__} for {lc.lasf} on pixel ({x}, {y}): \n")
                    print(f"\t{e}\n") 
    
        end = time.time()
        
        print(f'Metrics computed, time elapsed: {end - start}\n')
        
        ### STEP 3: SAVE VOXEL & GRID METRICS
        
        # Save outputs as pickles
        # "Can't open a pickle you don't know" - there can be malicious pickles, be wary.
        with open(f'{od_metrics}/{projstr}_{xysize}mgrid_covermetrics.obj', 'wb') as of:
            pickle.dump(lc.cover_dict, of, protocol=pickle.HIGHEST_PROTOCOL)
    
        with open(f'{od_metrics}/{projstr}_{xysize}mgrid_percmetrics.obj', 'wb') as of:
            pickle.dump(lc.perc_dict, of, protocol=pickle.HIGHEST_PROTOCOL)
    
        with open(f'{od_metrics}/{projstr}_{xysize}mgrid_heights.obj', 'wb') as of:
            pickle.dump(lc.height_dict, of, protocol=pickle.HIGHEST_PROTOCOL)
            
        with open(f'{od_metrics}/{projstr}_{xysize}mgrid_complexitymetrics.obj', 'wb') as of:
            pickle.dump(lc.complexity_dict, of, protocol=pickle.HIGHEST_PROTOCOL)
        
        # DONE
        endproj = time.time()
        projtime = endproj - startproj
        
        print(f'{lc.las.header.point_count} points gridded into {numcells} {xysize}m pixels in {projtime} seconds!\n')
        
        print(f'\tDone with processing {projstr}.\n')
        
    #===================================================================
    # Quick Script for Benchmarking time
    
    # start = time.time()
    
    # # set the cell indices to loop over in parallel
    # indices = lc.grid_dict['idx_cells'][0:1000]
    
    # t = []
    
    # ## Use concurrent futures to compute cover over each cell in parallel
    # with concurrent.futures.ProcessPoolExecutor(max_workers=None) as executor:
    #     for cphc, x, y in zip(executor.map(calccover_parallel, indices),
    #                              lc.grid_dict['x_cells'],
    #                              lc.grid_dict['y_cells']):
    
    #         try:
                
    #             # Stick the cover, perc, and heights inside the metrics dictionary
    #             # with x and y location as tuple keys
    #             lc.cover_dict[(x, y)] = cphc[0]
    #             lc.perc_dict[(x, y)] = cphc[1]
    #             lc.height_dict[(x, y)] = np.round(cphc[2], decimals=3)
    #             lc.complexity_dict[(x, y)] = cphc[3]
    
    #         except Exception as e:
    
    #             print(f"Saving metrics error - {e.__class__} for {lc.lasf} on pixel ({x}, {y}): \n")
    #             print(f"\t{e}\n") 
                
    # print(f'\t1000 {xysize} m cells takes {time.time() - start} s')
    # # Note: 1000 0.5 m cells takes 43.7 seconds
    # So you can't run this over an entire 1 km tile (~4 million 0.5 m pixels)
    # that would take about 2 days per tile.
    # running at 1 m pixels takes about half the time (# Note: 1000 1 m cells takes 25.5 seconds)
    # you have to run it over a smaller section, 
    # or downsize the resolution (that said, 1000 5 m cells takes 24.8 s, so there's not much gain)
    # either that, or make more of an effort to make the code efficient.
    #===================================================================
    #===================================================================
#tlas_voxelize(las_inputs,las_directory,output_directory,max_heigh######t,xysize,units="m"):

las_inputs=['/home/lefsky/time_trials/tile_66_136_sub_04.laz']
lc,index,struct=tlas_voxelize(las_inputs,"/home/lefsky/time_trials/","/home/lefsky/time_trials/",25,1,units="f")
    # Import Dependencies

