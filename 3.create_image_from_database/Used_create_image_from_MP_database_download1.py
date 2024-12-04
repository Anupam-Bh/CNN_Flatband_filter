# -*- coding: utf-8 -*-
"""
Created on Thu Jan 25 15:18:38 2024
Create image from band structure data of database 
@author: Anupam
"""
primarypath= 'D:/MatProj/electronic_str/pics/'
savepath='pics/'
my_dpi=120
#%%
import pickle
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from PIL import Image
import os
import gc
#import tracemalloc
#%%  
##### Not required CELL
#make one image joining multiple
def join_multiple_images(images):
    widths, heights = zip(*(i.size for i in images))
    total_width = sum(widths)
    max_height = max(heights)
    new_im = Image.new('L', (total_width, max_height))
    x_offset = 0
    for im in images:
      new_im.paste(im, (x_offset,0))
      x_offset += im.size[0]
    return new_im

#delete files
def delete_images(filename):
    for f in filename:
        fname = f.rstrip()
        if os.path.isfile(primarypath+fname): # this makes the code more robust
            os.remove(primarypath+fname) 
            
#%%
def extract_BS_image(bands_dict,matid,primarypath):
    ########Rescales each k-branch into unifrom width of 96 pixels
    eigvals_down=[]
    eigvals_up=[]
    arr_up=[]
    arr_down=[]
    #mat_structure=bands_dict.structure
    spinpol=bands_dict.is_spin_polarized
    branches=[list(x.values())[2] for x in bands_dict.branches]
    nodes1=[list(x.values())[0] for x in bands_dict.branches]
    nodes2=[list(x.values())[1]+1 for x in bands_dict.branches]  ### not to leave any gap
    #print(matid,nodes1,nodes2,branches)
    if len(branches)==1:
        return
    ef=bands_dict.efermi
    #print('Fermi energy = ',ef, '\n Material formula =',mat_structure.alphabetical_formula)
    eigvals_up=list(bands_dict.bands.values())[0]
    if spinpol is True:
        eigvals_down=list(bands_dict.bands.values())[1]
        #print(eigvals_up.shape, eigvals_down.shape)
    filename=[]
    splitpath=splitpath = [i.split('-') for i in branches]
    plotpath=[]
    path_idx=[]
    for num,paths in enumerate(splitpath):
        if paths[0]!=paths[1]:
            plotpath.append(paths)
            path_idx.append(num)
    fi,axxs =plt.subplots(nrows=1, ncols=len(plotpath), figsize=((1.0001*96*len(plotpath))/my_dpi, (4*96)/my_dpi), dpi=my_dpi)  
    for num,pathid in enumerate(path_idx):
        arr_up=eigvals_up[:,nodes1[pathid]:nodes2[pathid]]
        if spinpol==True:
            arr_down=eigvals_down[:,nodes1[pathid]:nodes2[pathid]]
        for i in range(arr_up.shape[0]):
            axxs[num].plot(arr_up[i,:],color='black',linestyle='-',linewidth=2.5) 
            if spinpol == True:
                axxs[num].plot(arr_down[i,:],color='black',linestyle='-',linewidth=2.5) 
            axxs[num].set_ylim([-1+ef,1+ef])
            axxs[num].axis('off')
            axxs[num].margins(0,0)
    fid = matid+'.png'
    fi.tight_layout()
    plt.subplots_adjust(left = 0, right =1.0, top = 1.0, bottom = 0, wspace=0, hspace=0)
    fi.savefig(savepath+fid, bbox_inches='tight', pad_inches = 0)
    plt.clf()
    plt.close('all')    ### crucial
    
    
#%%
def unpack_each_file_and_plot(startid,endid,first_id,last_id,Masterlist):
    #gc.collect()
    with open(primarypath+'../'+'BS_MP_'+str(startid)+'_'+str(endid),'rb') as op:
        [bs,ids]= pickle.load(op)
    print(len(bs))
    print(len(ids))
    #gc.collect()
    for mat_counter,bands_dict in enumerate(bs):
        if int(ids[mat_counter].split('-')[1]) >= first_id and int(ids[mat_counter].split('-')[1]) <= last_id and (ids[mat_counter] in Masterlist):
            extract_BS_image(bands_dict,ids[mat_counter],primarypath)
            #print(ids[mat_counter], bands_dict.branches)
            if mat_counter % 10 ==0:
                print(ids[mat_counter])
    del bs, ids, bands_dict
#%%   Create images for required Material IDs which are in the Master list
##Masterlist
with open('list_of_materials_with_BS_DOS','rb') as ff:
    Masterlist= pickle.load(ff)
##
first_id= 772755
last_id=  890000
list_BS_files=glob.glob(primarypath+'../'+'BS_MP_*')
List_file_required=[]

for i in list_BS_files:
    fname=i.split('\\')[1]
    startid = fname.split('_')[2]
    endid = fname.split('_')[3]
    if float(startid) <= last_id and float(endid) >= first_id:
        List_file_required.append(i)
        print(List_file_required[-1])
        unpack_each_file_and_plot(startid,endid,first_id,last_id,Masterlist)
        
        
