#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 18:21:33 2024

@author: ubuntu
"""

import numpy as np
from itertools import chain
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import sys

# # Create list of images
file_list=['nodup_100_22_density_lt1','nodup_100_22_density','nodup_100_22_cover',
 'nodup_100_22_pdtm','nodup_100_22_pdsm','nodup_100_22_chm','nodup_512_11_density_lt1',
 'nodup_512_11_density','nodup_512_11_cover','nodup_512_11_pdtm','nodup_512_11_pdsm','nodup_512_11_chm',
 'nodup_669_8_density_lt1','nodup_669_8_density','nodup_669_8_cover','nodup_669_8_pdtm','nodup_669_8_pdsm','nodup_669_8_chm',
 'nodup_985_4_density_lt1','nodup_985_4_density','nodup_985_4_cover','nodup_985_4_pdtm','nodup_985_4_pdsm','nodup_985_4_chm',
 'sub_01_chm','sub_01_cover','sub_01_density_lt1','sub_01_density','sub_01_pdsm','sub_01_pdtm']

title=['Density_lt1','Density','Cover (%)','Maximum Canopy Elevation (ft)','Terrain Elevation (ft)','Canopy Height (ft)',
       'Density_lt1','Density','Cover (%)','Terrain Elevation (ft)','Terrain Elevation (ft)','Canopy Height (ft)',
       'Density_lt1','Density','Cover (%)','Terrain Elevation (ft)','Terrain Elevation (ft)','Canopy Height (ft)',
       'Density_lt1','Density','Cover (%)','Terrain Elevation (ft)','Terrain Elevation (ft)','Canopy Height (ft)',
       'Density_lt1','Density','Cover (%)','Terrain Elevation (ft)','Terrain Elevation (ft)','Canopy Height (ft)']



#print(len(file_list))
# print(len(file_list))
reference_bands=['sub_01_chm','sub_01_cover','sub_01_density_lt1','sub_01_density','sub_01_pdsm','sub_01_pdtm']
 
import rasterio

from collections.abc import Iterable
    
def flatten(coll):
    for i in coll:
            if isinstance(i, Iterable) and not isinstance(i, str):
                for subc in flatten(i):
                    yield subc
            else:
                yield i

def unpack_dict_var(dict_array,varname):
 #   print(type(dict_array),varname)
    all=[]
 #   print(type(dict_array[0]))

    for x in dict_array:
#        print(x)
        all.append(x[varname])
        
 #   all=[print(x) for x in dict_array]
#    print("upd_sub",all,np.array(all))
    return(all)

def extract_fileinfo(filenames):
    out=[]
    for ix,f in enumerate(filenames):
        fs=f.split("_")
        if fs[-1] == "lt1":
            fs=[fs[0:len(fs)-2],fs[-2]+"_lt1"]
            fs=list(flatten(fs))
        if len(fs)==3:
            fs=[fs[0],fs[1],str(100),fs[2]]
#        print(fs)
        fs=list(flatten(fs))

        tmp={"filename":f,f"ilename_id":'_'.join(fs),"shot_density":fs[2],"var":fs[3],"ix":ix}     
        out.append(tmp)
    return(out)
    
ef=extract_fileinfo(file_list)
#
#print(ef[0])

#==========================================================================
#
#
#  You edited the first of the next two accidentaly so you'll need to 
#  reconstruct the algorithm from both
#
#  Good luck!
#
#==========================================================================

# def match_and_substract(efinfo,in_image,titles):
    
#     all_filenames=unpack_dict_var(efinfo,'filename')    
#     all_shot_density=unpack_dict_var(efinfo,'shot_density')
#     all_var=unpack_dict_var(efinfo,'var')
#     all_index=unpack_dict_var(efinfo,'ix')

#     baseline_filenames=all_filenames[all_shot_density == str(100)]
#     baseline_vars=all_var[all_shot_density == str(22)]
#     baseline_index=all_index[all_shot_density == str(22)]
    
    
#     infile_raster = rasterio.open(in_image)
#     mask=infile_raster.read_masks()
#     mask=mask[0,:,:].squeeze()
#     mask[:]=1
#     print(mask.shape)
#     meta = infile_raster.meta.copy()
#     out_meta=infile_raster.meta.copy()
#     out_meta.update({"count":out_meta["count"]*2})

        
#     with rasterio.open(in_image.replace(".tif","_wdiff.tif"), 'w', **out_meta) as dest:
# #        print(type(dest))
# #        print(len(list(range(1,meta["count"]+1))))
# #        print(len(efinfo))
#         print(all_index)
#         for ix in all_index:
#             print("ix",ix)
#             r=infile_raster.read(ix)
#             dest.write(infile_raster.read(ix), ix)
#             dest.set_band_description(ix, file_list[ix-1])
#             af=efinfo[ix-1]
#             match=baseline_filenames[af['var'] == baseline_vars]
#             match_ix=baseline_index[af['var'] == baseline_vars]
#             tmp1=infile_raster.read(ix).squeeze()
#             tmp2=infile_raster.read(match_ix+1).squeeze()
#             out_calc=tmp1-tmp2
#             mask=mask*(tmp1>-9999)*(tmp2>-9999)
#             dest.write(out_calc, ix+len(efinfo))
#             dest.set_band_description(ix+len(efinfo), 'diff_'+file_list[ix-1])
#             rtmp=r[r > -9999]
#             rtmp2=out_calc[np.logical_and(tmp1 > -9999,tmp2 > -9999) ]
# #            print(rtmp2)
# #            rtmp2=r[r >= -9999]
#             print(ix,file_list[ix-1],np.mean(rtmp),min(rtmp),max(rtmp),np.std(rtmp),"diff_"+file_list[ix-1],np.mean(rtmp2),min(rtmp2),max(rtmp2),np.std(rtmp2))
#             plt.hist(rtmp2,bins=50,weights=np.ones(len(rtmp2)) / len(rtmp2))
#             plt.title("Difference in 8ppm and 22ppm estimates of "+titles[ix-1],fontsize=10)
#             plt.xlabel('Differrence',fontsize=10)
#             plt.ylabel('Percent of observations',fontsize=10)   
#             plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
#             plt.show()
#             print(ix)
            
# #            input("Press Enter to continue...")
                    
#         print(np.min(mask),np.max(mask))
#         dest.write_mask(mask*255)

# #   for ix,af in enumerate(efinfo):

     
#     print(mask.shape)
     
#     return(baseline_filenames,baseline_vars,af)

def return_struct_list_wcriteria(instruct,variable,desired_shotdens):
#    for i in instruct:
#        print(i['shot_density'])
    shot_density=np.array(unpack_dict_var(instruct,'shot_density'))
    print("SDSDSD",shot_density)
#    print(variable)
#    print(desired_shotdens)
    if desired_shotdens == None:
        good_ix=list(range(0,len(shot_density)))
        print("GOOD_IDX1",good_ix)

    else: 
        good_ix=list(np.where(shot_density == desired_shotdens))
        print(shot_density)
        print("GOOD_IDX2",good_ix)
    out=[]
#    print(">>>3",good_ix[0])
#    print(">>>2",len(good_ix[0]))
#    print(">>>1",len(instruct))
    print(good_ix,len(good_ix))
    for x in good_ix[0]:
        print(">>>",x,variable,instruct[x])
        out.append(instruct[x][variable])
    return(out)

def match_and_substract(efinfo,in_image,titles,verbose=False):#efinfo
    
 #   print(efinfo[0].keys())
##    comparison_shot_density="8"
#    all_shot_density=return_struct_list_wcriteria(efinfo,"shot_density",comparison_shot_density)
#    all_filenames=return_struct_list_wcriteria(efinfo,'filename',comparison_shot_density) 
#    all_vars=return_struct_list_wcriteria(efinfo,'var',comparison_shot_density)
#    all_index=list(range(0,len(all_vars)))
#    all_index2=return_struct_list_wcriteria(efinfo,'ix',comparison_shot_density)#

    all_shot_density=unpack_dict_var(efinfo,'shot_density')
    all_filenames=unpack_dict_var(efinfo,'filename')
    all_vars=unpack_dict_var(efinfo,'var')
    all_index=list(range(0,len(all_vars)))
    all_index2=unpack_dict_var(efinfo,'ix')
    #print("alix",all_index)
    #print(all_index2)
#    print("LBA: ",len(all_index),len(all_index),len(all_vars),len(all_filenames),len(all_shot_density))

#    print(all_shot_density == str(22))
#    print(type(all_shot_density == str(22)),type(baseline_filenames))
    
    baseline_target_density="22"
    baseline_filenames=return_struct_list_wcriteria(efinfo,'filename',baseline_target_density) 
    baseline_vars=return_struct_list_wcriteria(efinfo,'var',baseline_target_density)#
    baseline_shot_density=return_struct_list_wcriteria(efinfo,'shot_density',baseline_target_density)
    baseline_index=list(range(0,len(baseline_vars)))
    baseline_index2=return_struct_list_wcriteria(efinfo,'ix',baseline_target_density)

    #print("LBA: ",len(baseline_index),len(all_index),len(baseline_vars),len(baseline_filenames),len(baseline_shot_density))
    
                                
    infile_raster = rasterio.open(in_image)
    mask=infile_raster.read_masks()
    mask=mask[0,:,:].squeeze()
    mask[:]=1
#    print(mask.shape)
    meta = infile_raster.meta.copy()
    out_meta=infile_raster.meta.copy()
    out_meta.update({"count":out_meta["count"]*2})

    if verbose:
        print("ai",all_index)
        print("bi",baseline_index)
    with open("/home/ubuntu/density_comparison_v2/difference_output_"+baseline_target_density+".txt", 'w') as f:
        header="shot_density,variable,all_filenames,all_mean,all_std,all_min,all_max,base_filenames,base_mean,base_std,base_min,base_max,diff_mean,diff_std,diff_min,diff_max"
        f.write(header+"\n")
        with rasterio.open(in_image.replace(".tif","_wdiff.tif"), 'w', **out_meta) as dest:
             for ix in all_index:
    #            ix=ix+1match_ixz
                
                if verbose:
    
                    print("==========================================================================")
                    print("==========================================================================")
                    print("==========================================================================")
                
                    print("IX",ix)
    
                r=infile_raster.read(ix+1)
                dest.write(r, ix+1)
                
                af=efinfo[ix]
    #            print("aff",af)
    #            print(baseline_var)
       #            print("afb",af.shape,baseline.shape)
                #match=baseline_vars[af['var'] == baseline_vars]
                if verbose:
                    print('len',len(baseline_vars),len(baseline_index),len(af))
                    print('type',type(baseline_vars),type(baseline_index),type(af))
    
    #            print("-----------------------------------------------------",baseline_index)
    #            print(af['var'],np.array(baseline_vars))
    #            for mix in range(baseline_index):
                    #print("compare var",af['var'])
                    #print("bv",baseline_vars,len(baseline_vars),len(baseline_index))
                    #print("av",all_vars,len(all_vars),len(all_index))
                    #print()
    
    
    
    
    
                match_ix=[]
                match=[]
                for b in range(len(baseline_vars)):
                    #print("CHECK",baseline_vars[b],af['var'],baseline_vars[b]==af['var'])
                    if (baseline_vars[b]==af['var']):
                        match_ix.append(b)
                        match.append(baseline_filenames[b])
                #print("MATCH",match)
                try:
                     tmp=len(match_ix)
                except:
                     print("MATCH_ISNT_LISt")
                else:            
                    #print(tmp,len(match_ix))
                    match_ix=match_ix[0]
    
    
    
    
                if verbose:
                    print("111111111111111111111111111111111111111111111111111111111111111111111111111")
                    print("TRYTRI",match_ix,baseline_vars[match_ix],af['var'])
                    print(len(all_filenames),len(baseline_filenames))
                    print("--------------------")
                    print("MATCHIX",match_ix)
                    print("KATCH",match)
                    print(af['var'])
                    print("IX_MATCH",ix,match_ix)
                    print("22222222222222222222222222222222222222222222222222222222222222222222222222222")
                    print("<<<<<",ix,match_ix)
                    print(baseline_filenames[match_ix],all_filenames[ix])
                    print(">>>>>",ix,match_ix,baseline_vars[match_ix],all_vars[ix])
                
                    print([af['var']])
                    print('MATCH:2',match,match_ix,ix,baseline_filenames[match_ix],all_filenames[ix])  
                    print(baseline_index[match_ix],all_index[ix])
    
                all_index2=list(all_index2)
    
                #print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")
                #print(all_index2[ix]+1,all_index2)
                #print(len(all_index2))
                #print(type(all_index2))
                #print(type(ix))
                #print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")
                #print(baseline_index2[ix]+1,baseline_index2)
                #print(len(baseline_index2))
                #print(type(baseline_index2))
                #print(type(ix))
                #print(ix)
    #                  [ix]+1),type(all_index2),type(ix))            ##print(len(all_index2[ix]+1))#,
    
      #          print(">>>>>",    all_index2[int(match_ix)+1],
#                print(ix,len(all_index2),match_ix,len(baseline_index2))
                tmp1=infile_raster.read(all_index2[int(ix)]+1).squeeze()
                tmp2=infile_raster.read(baseline_index2[int(match_ix)]+1).squeeze()
                
                out_calc=tmp2-tmp1
                if verbose:
                    print("33333333333333333333333333333333333333333333333333333333333333333333333333333")  
                    print(baseline_vars)
                    print("band_stats",np.mean(tmp1),np.mean(tmp2))
    #            mask=mask*(tmp1>-9999)*(tmp2>-9999)
                dest.write(out_calc, ix+len(efinfo))
                dest.set_band_description(ix+len(efinfo), 'diff_'+file_list[ix-1])
                if (af['var'] !='cover'):
                    rtmp=r[r >1]
                else:   
                    rtmp=r[r>0]
                if verbose:
                    print("444444444444444444444444444444444444444444444444444444444444444444444444444444")
                    print("SHAPE",tmp1.shape)
                is_good=np.logical_and(tmp1 > -9999,tmp2 > -9999)
                tmp1=tmp1[is_good]
                tmp2=tmp2[is_good]
                out_calc=tmp2-tmp1
                rtmp2=out_calc
#                print(len(tmp1),len(tmp2),len(rtmp2))

                np.set_printoptions(linewidth=200)
                
                tmp=all_filenames[ix].split("_") 
#                      print(tmp)
                if len(tmp) == 5:
                  var_out=tmp[3]+"_"+tmp[4]
                  dns_out=tmp[2]
                if len(tmp)==4:
                  var_out=tmp[3]
                  dns_out=tmp[2]
                if len(tmp)==3:
                  var_out=tmp[2]
                  dns_out="100"
       
#                print(var_out)
#                print(dns_out)
                          
 #               print(len(var_out,len(dns_out),len(all_filenames))
 #               print(lvar_out),len(dns_out),len(all_filenames))
                output1=[all_filenames[ix],var_out,dns_out,np.mean(tmp2),np.std(tmp2),min(tmp2),max(tmp2),
                      baseline_filenames[match_ix],np.mean(tmp1),np.std(tmp1),min(tmp1),max(tmp1),
                      np.mean(rtmp2),np.std(rtmp2),min(rtmp2),max(rtmp2)]
                output=[]
                for o in output1:
                    output.append(str(o))
#                output=','.join(str(v) for v in output)
                print(str(",").join(output))
                f.write(",".join(output)+"\n")
    #            print("aff",match,af['filename'])
                #print("1     ","Diff   ","mean(rtmp2),np.std(rtmp2),min(rtmp2),max(rtmp2))
                #print("2     ",baseline_filenames[match_ix],np.mean(tmp1),np.std(tmp1),min(tmp1),max(tmp2))
                #print("3     ",all_filenames[ix ],np.mean(tmp2),np.std(tmp2),min(tmp2),max(tmp2))
     #           #print ("diff_"+file_list[ix-1],np.mean(rtmp2),min(rtmp2),max(rtmp2),np.std(rtmp2))
                xrange=(np.percentile(rtmp2,0.5),np.percentile(rtmp2,99.5))
                plt.hist(rtmp2,bins=40,range=xrange,density=True)#,weights=np.ones(len(rtmp2)) / len(rtmp2))
                plt.title("8 ppm and 22 ppm estimates of "+titles[ix-1],fontsize=10)
                plt.xlabel('Difference',fontsize=10)          
                plt.ylabel('Percent of observations',fontsize=10)   
                plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
                if dns_out == "8":
                    plt.savefig('/home/ubuntu/density_comparison_v2/diffig_'+var_out+".png")
                    print("savefig")
#                plt.show()
#            #print(ix)
            
#            input("Press Enter to continue...")
                    
#        print(np.min(mask),np.max(mask))
             dest.write_mask(mask*255)

#   for ix,af in enumerate(efinfo):

     
    print(mask.shape)
     
    return(baseline_filenames,baseline_vars,af)

print(ef[0])

#mas=match_and_substract_8_22(ef,"/home/ubuntu/density_comparison_v2/multistack2.tif",title)
mas=match_and_substract(ef,"/home/ubuntu/density_comparison_v2/multistack2.tif",title)


