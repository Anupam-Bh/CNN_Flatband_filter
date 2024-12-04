import json
#import sys
import os
#import pymatgen
from pymatgen.core import Structure #, Molecule, Lattice
from matminer.featurizers.site import CrystalNNFingerprint
from matminer.featurizers.structure import SiteStatsFingerprint
import numpy as np
#import matplotlib.pyplot as plt
#import chemparse
import glob
import pickle
#### Self defined class
os.chdir('C:/Users/Anupam/Desktop/MP_flatband/flatness_code/')
from flatband_materials_from_MP import flatband_materials_from_MP
path = 'D:/MatProj/electronic_str/'


#%% define fingerprints
ssf = SiteStatsFingerprint(CrystalNNFingerprint.from_preset('ops', distance_cutoffs=None, x_diff_weight=0),stats=('mean','std_dev','maximum','minimum'))
#ssf = SiteStatsFingerprint(CrystalNNFingerprint.from_preset('cn'),stats=('mean'))

#%%  Read the required DOS files from MP database and append
def read_DOS_files(path,first_id,last_id):
    DOS_all= []
    id_all = []
    list_DOS_files=glob.glob(path+'DOS_MP_*')
    List_DOS_file_required=[]
    for i in list_DOS_files:
        fname=i.split('\\')[1]
        startid=fname.split('_')[2]
        endid=fname.split('_')[3]
        if float(startid) <= last_id and float(endid) >= first_id:
            List_DOS_file_required.append(i)
            with open(i,'rb') as oop:
                [DOS,ids]=pickle.load(oop)
            print('Reading DOS file:'+ str(i))
            DOS_all.extend(DOS)
            id_all.extend(ids)
    return DOS_all, id_all

#%%      Truncating the structure from its basis: 
    ################    only one type of elements remains based on DOS or min stoichiometry in case DOS NA
def trunacte_structure(matid,dos_file,bb):
    aa=dos_file.structure.as_dict()
    #print(type(aa['sites']))
    sites=[]
    for n in aa['sites']:
        #print(n)
        if n['label']==bb:
            #print(n)
            sites.append(n)
        #else:
         #   aa['sites'].remove(n)
    #print(sites)
    aa['sites']=sites
    #print(aa)
    #############################################
    structure = Structure.from_dict(aa)
    #print(structure)
    return structure
#%% get fingerprint dictionary
def get_fingerprint(path, first_id, last_id, Mat_flat):
    sublattice_wrong=0
    D={}
    [DOS_all,id_all] = read_DOS_files(path,first_id,last_id)
    for obj in Mat_flat:
        matid=obj.material_id
        if (matid in id_all) and DOS_all[id_all.index(matid)].get_element_dos()  and obj.flatness: # when (Flat) and element projected DOS available
            
            ### Find species which has largest DOS
            dos_file = DOS_all[id_all.index(matid)]  ### when atom projected DOS is available
            elems = dos_file.structure.elements
            elem_DOS = dos_file.get_element_dos()
            ldos_max={}
            ldos_sum={}
            ef=dos_file.efermi
            eng=dos_file.energies
            for i in elems:
                denst= elem_DOS[i].as_dict()['densities']
                for j in denst.keys():   #for spin polarized case
                    if j=='1':
                        sum_density=np.array(denst[j])
                    else:
                        sum_density=sum_density+np.array(denst[j])
                dos_full = np.stack((eng-ef,sum_density),axis=1)
                dos_seg=  obj.flat_segments_preferred()      ####find which horizontal segment of DOS is chosen
                range_dos=[1.5-(dos_seg*0.5),1-(dos_seg*0.5)]
                b = dos_full[dos_full[:, 0] <= range_dos[0]]
                b = b[b[:, 0] >= range_dos[1]]
                #print(b)
                ldos_max[i]=max(b[:,1])
                ldos_sum[i]=np.trapz(b[:,1],b[:,0])
            bb=max(ldos_max,key=ldos_max.get)    ### sublattice found from maximum peak height in bandwidth
            bb_sumdensity=max(ldos_sum,key=ldos_sum.get)    ### sublattice identified from heighest integrated PDOS in bandwidth
            #print(dos_max[i])
            ldos_sum=dict(sorted(ldos_sum.items(), key=lambda item: item[1],reverse=True))
            truncated_structure =trunacte_structure(matid,dos_file,str(bb)) 
            
            
            
            ### Find if predicted sublattices are same
            if bb!=bb_sumdensity:
                #print(matid+' '+str(elems))
                #print('sublattice prediction wrong')
                sublattice_wrong=sublattice_wrong+1
            
            
            
            ### Calculate sublattice scores       
            score_prim=0
            score_sec=0
            score_ter=0
            score_prim=list(ldos_sum.values())[0]/sum(ldos_sum.values())
            #print(score_prim)
            if len(ldos_sum)>1:
                score_sec=list(ldos_sum.values())[1]/sum(ldos_sum.values())
            if len(ldos_sum)>2:
                score_ter=list(ldos_sum.values())[2]/sum(ldos_sum.values())
                #print(score_sec)
            score_max=max(ldos_max.values())/sum(ldos_max.values())
            score_sum=max(ldos_sum.values())/sum(ldos_sum.values())
            print(matid+' '+str(elems)+ '  Flatband sublattice:'+str(bb)+'  Primay score:'+str(score_prim)+'  '+str(list(ldos_sum.values())[0])+'  '+ str(sum(ldos_sum.values())))
            
            ### Create dictionary to store the information of sublattice, fingerprint and scores
            #print(str(bb),score_prim, dos_file.structure.as_dict(), truncated_structure, dos_file.structure.formula, elems)
            try:
                vector=ssf.featurize(truncated_structure)
                D[matid]={}
                D[matid]['structure_fingerprint']=vector
                D[matid]['sublattice_element']= str(bb)
                D[matid]['primary_sublattice_score']=score_prim
                if len(ldos_sum)>1:
                    D[matid]['secondary_sublattice_score']=score_sec
                if len(ldos_sum)>2:
                    D[matid]['tertiary_sublattice_score']=score_ter
                D[matid]['full_structure'] = dos_file.structure.as_dict()
                D[matid]['truncated_structure'] = truncated_structure.as_dict()
                D[matid]['formula'] = dos_file.structure.formula
                #D[matid]['emements'] = elems
                D[matid]['Fermi'] = ef
                
                del matid, vector, score_prim, score_sec, score_ter, score_max, score_sum
            except:
                print('Error:  '+matid+'     '+dos_file.structure.formula)
        elif obj.flatness:         #when no DOS or no atom projected DOS
            print("No DOS for %s" % matid)
        else:
            pass
            #print("No Flatband found for %s" % matid)
    
    
    with open('CrystalNN_fingerprint_with_structure_from_DOS_'+str(first_id)+'_'+str(last_id)+'.json', 'w', encoding='utf-8') as f:
            json.dump(D, f, ensure_ascii=False, indent=4)
    print(str(sublattice_wrong) +' sublattice prediction mismatch when used sum of density instead of max density')



#%%  Read flatness flatness prediction
list_prediction_files = glob.glob('4_predict/list_predicted_materials_*')
starts=[]
ends=[]
for items in list_prediction_files:
    starts.append(int(items.split('_')[-2]))
    ends.append(int(items.split('_')[-1]))
starts.sort()
ends.sort()
for i in range(0,len(starts)-34):   # from 1 to rest
    first_id=starts[i]
    last_id= ends[i]
    with open('4_predict/list_predicted_materials_'+str(first_id)+'_'+str(last_id),'rb') as op:
        Mat_flat= pickle.load(op)
    print('Mat_flat reading done')
    get_fingerprint(path, first_id, last_id,Mat_flat)
    










