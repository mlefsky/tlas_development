#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 23:43:23 2024

@author: ubuntu



las_directory="/home/ubuntu/density_comparison_v2/"
dummy=image_to_table(las_directory+"chm_stack.tif",10000,
                     las_directory+"chm_stack.csv",3.28)



las_directory="/home/ubuntu/density_comparison_v2/"
dummy=image_to_table(las_directory+"pdsm_stack.tif",10000,
                     las_directory+"pdsm_stack.csv",3.28)



las_directory="/home/ubuntu/density_comparison_v2/"
dummy=image_to_table(las_directory+"pdtm_stack.tif",10000,
                     las_directory+"pdtm_stack.csv",3.28)




las_directory="/home/ubuntu/density_comparison_v2/"
dummy=image_to_table(las_directory+"density_stack.tif",10000,
                     las_directory+"density_stack.csv",3.28)

las_directory="/home/ubuntu/density_comparison_v2/"
dummy=image_to_table(las_directory+"density_lt1_stack.tif",10000,
                     las_directory+"density_lt1_stack.csv",3.28)





las_directory="/home/ubuntu/density_comparison_v2/"
prefix_list=["nodup_100_22","nodup_512_11","nodup_985_4","nodup_669_8"]

for prefix in prefix_list:
    
    dummy=image_to_table(las_directory+prefix+"_stack.tif",10000,
                         las_directory+prefix+"_stack.tif",3.28)




las_directory="/home/ubuntu/density_comparison_v2/"
dummy=image_to_table(las_directory+"nodup_100_22_stack.tif",10000,
                     las_directory+"nodup_100_22_stack.csv",3.28)

las_directory="/home/ubuntu/density_comparison_v2/"
dummy=image_to_table(las_directory+"nodup_512_11_stack.tif",10000,
                     las_directory+"nodup_512_11_stack.csv",3.28)

las_directory="/home/ubuntu/density_comparison_v2/"
dummy=image_to_table(las_directory+"nodup_985_4_stack.tif",10000,
                     las_directory+"nodup_985_4_stack.csv,3.28)

las_directory="/home/ubuntu/density_comparison_v2/"
dummy=image_to_table(las_directory+"nodup_669_8_stack.tif",10000,
                     las_directory+"nodup_669_8_stack.csv",3.28)



