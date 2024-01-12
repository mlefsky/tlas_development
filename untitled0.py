#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 15:51:52 2024

@author: mlefsky
"""

#tmp=main()
plt.clf()
v1,v2,d=tmp
v1_maxz=plot_voxel_slice(v1['voxels'][250:550,200,:],v1,title='Transect',show_image=True)
v2_maxz=plot_voxel_slice(v2['voxels'][250:550:,200,:],v1,title='Transect',show_image=False)
plt.plot(np.flip(v2_maxz),color="r",linewidth=4,markersize=24)
plt.plot(np.flip(v1_maxz),color="b",linewidth=4)
plt.plot(np.flip(v1_maxz-v2_maxz))
print(np.sum(v1['voxels'][:,200,:]))
print(np.sum(v2['voxels'][:,200,:]))
