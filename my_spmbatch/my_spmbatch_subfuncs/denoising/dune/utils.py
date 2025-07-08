#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 14:36:39 2024

@author: dr. Peter Van Schuerbeek (UZZ Brussel)
"""

import os
import json
import math
import nibabel as nib
import numpy as np

from nilearn.image import mean_img, clean_img
from scipy.stats.mstats import zscore
from scipy.signal import resample

"""
---------------------------------------------------------------------------------------
"""     
    
def get_brain_and_non_brain_masks(datpath,gmfile,wmfile,csffile,maskfile):
    
    """
    Make brain (GM+WM) and no-brain masks 
    """
    
    gmim = nib.load(gmfile)
    wmim = nib.load(wmfile)
    csfim = nib.load(csffile)
    
    gm_data = gmim.get_fdata()
    gm_data = np.nan_to_num(gm_data,nan=0.0)
    
    wm_data = wmim.get_fdata()
    wm_data = np.nan_to_num(wm_data,nan=0.0)
    
    csf_data = csfim.get_fdata()
    csf_data = np.nan_to_num(csf_data,nan=0.0)
    
    maskim = nib.load(maskfile)  
    gmwmdata = maskim.get_fdata()
    gmwmdata[gm_data+wm_data<0.8]=0
    gmwmdata[gmwmdata>0.0]=1
                
    maskim = nib.load(maskfile)
    nbdata = maskim.get_fdata()
    nbdata[csf_data<0.8]=0
    nbdata[gm_data+wm_data>0.1]=0
    nbdata[nbdata>0.0]=1
    
    maskim = nib.load(maskfile)
    braindata = maskim.get_fdata()
    braindata[gm_data+wm_data+csf_data<0.1]=0
    braindata[braindata>0.0]=1
            
    print('make masks done')
    
    return gmwmdata, nbdata, braindata

"""
---------------------------------------------------------------------------------------
""" 

def get_func_data(datpath,func,jsonfile,nechoes):
    
    """
    Load brain (GM+WM) and no-brain data 
    """
    
    for e in list(range(nechoes)):
        
        if nechoes>1:
            splitf = func.split('_echo-1')
            tfunc = datpath+splitf[0]+'_echo-'+str(e+1)+splitf[1]
        else:
            tfunc = datpath+func
        
        tfmridat = nib.load(tfunc)     
        tfdat = tfmridat.get_fdata()  
        
        tfdat = np.transpose(np.nan_to_num(tfdat,nan=0.0),axes=(3,0,1,2))
        
        if e==0:
            fdat = np.reshape(tfdat,tfdat.shape+(1,))
        else:
            fdat = np.append(fdat,np.reshape(tfdat,tfdat.shape+(1,)),axis=-1)
            
    print('get functional data done')
    
    return fdat

"""
---------------------------------------------------------------------------------------
""" 

def get_train_data(fdat,gmwmdata,nbdata,noise_file,jsonfile):
    
    """
    Make the data for training the model 
    """
    
    fdat = np.transpose(fdat,axes=(4,0,1,2,3))
    
    #Brain data to train the model
    gmfdat = fdat[:,:,gmwmdata>0.0].T
    
    nvoxel_gmfdat = gmfdat.shape[0]
    nvoxel_train = min(nvoxel_gmfdat,50000)
    
    trainind_gmfdat = np.random.permutation(int(nvoxel_gmfdat))[:int(nvoxel_train)]
    train_gm = gmfdat[trainind_gmfdat,:,:]
    
    #No brain data to train the model
    nbfdat = fdat[:,:,nbdata>0.0].T
    
    nvoxel_nbdat = nbfdat.shape[0]
    if nvoxel_nbdat<nvoxel_train:
        nbfdat = np.repeat(nbfdat, math.ceil(nvoxel_train/nvoxel_nbdat), axis=0)
        nvoxel_nbdat = nbfdat.shape[0]
    
    trainind_nbfdat = np.random.permutation(int(nvoxel_nbdat))[:int(nvoxel_train)]
    train_nb = nbfdat[trainind_nbfdat,:,:]
    
    #Noise data (e.g. motion regressors)
    params = np.genfromtxt(noise_file)
    noisedat = params.transpose()
    noisedat = np.reshape(noisedat,noisedat.shape+(1,))
    
    nvoxel_noisedat = noisedat.shape[0]
    if nvoxel_noisedat<nvoxel_train:
        noisedat = np.repeat(noisedat, math.ceil(nvoxel_train/nvoxel_noisedat), axis=0)
        nvoxel_noisedat = noisedat.shape[0]
        
    trainind_noisedat = np.random.permutation(int(nvoxel_noisedat))[:int(nvoxel_train)]    
    train_noise = noisedat[trainind_noisedat,:,:]
    
    #Standardize input
    train_noise = train_noise - np.mean(train_noise,axis=1,keepdims=True)
    train_noise = np.nan_to_num(train_noise / np.std(train_noise,axis=1,keepdims=True),nan=0.0)
    
    train_gm = train_gm - np.mean(train_gm,axis=1,keepdims=True)
    train_gm = np.nan_to_num(train_gm / np.std(train_gm,axis=1,keepdims=True),nan=0.0)
    
    train_nb = train_nb - np.mean(train_nb,axis=1,keepdims=True)
    train_nb = np.nan_to_num(train_nb / np.std(train_nb,axis=1,keepdims=True),nan=0.0)

    train_gm = np.transpose(train_gm,axes=(0,2,1))
    train_nb = np.transpose(train_nb,axes=(0,2,1))
    train_noise = np.transpose(train_noise,axes=(0,2,1))
   
    print('get training data done')
    
    return train_gm, train_nb, train_noise
    
"""
---------------------------------------------------------------------------------------
""" 

def get_apply_func_data(fdat,maskdata):
    
    """
    Prepare the data for the model 
    """
    
    fdat = np.transpose(fdat,axes=(4,0,1,2,3))
    
    func_data = fdat[:,:,maskdata>0.0].T
    
    func_data = func_data - np.mean(func_data,axis=1,keepdims=True)
    func_data = np.nan_to_num(func_data / np.std(func_data,axis=1,keepdims=True),nan=0.0)

    func_data = np.transpose(func_data,axes=(0,2,1))
    
    return func_data

"""
---------------------------------------------------------------------------------------
""" 

def save_denoised_data(dfdat,dnfdat,fdat,maskdata,func,datpath,isasl):
    
    """
    Save the denoised data
    """
    
    fdat = np.transpose(fdat,axes=(1,2,3,0,4))
    
    meanfdat = np.mean(fdat[maskdata>0.0,:,0],axis=1,keepdims=True)
    stdfdat = np.std(fdat[maskdata>0.0,:,0],axis=1,keepdims=True)
    
    meandfdat = np.mean(dfdat,axis=1,keepdims=True)
    stddfdat = np.std(dfdat,axis=1,keepdims=True)
    
    dfdat = np.nan_to_num(dfdat - meandfdat,nan=0.0)
    
    dfdat = np.nan_to_num(dfdat * stdfdat,nan=0.0)
    dfdat = np.nan_to_num(dfdat + meanfdat,nan=0.0)
    
    meandnfdat = np.mean(dnfdat,axis=1,keepdims=True)
    stddnfdat = np.std(dnfdat,axis=1,keepdims=True)
    
    dnfdat = np.nan_to_num(dnfdat - meandnfdat,nan=0.0)
    
    dnfdat = np.nan_to_num(dnfdat * stdfdat,nan=0.0)
    dnfdat = np.nan_to_num(dnfdat + meanfdat,nan=0.0)
    
    nfdat = np.zeros((fdat.shape[0],fdat.shape[1],fdat.shape[2],dfdat.shape[1]))
    nfdat[maskdata>0.0,:] = dfdat
    
    nnfdat = np.zeros((fdat.shape[0],fdat.shape[1],fdat.shape[2],dnfdat.shape[1]))
    nnfdat[maskdata>0.0,:] = dnfdat
    
    nfdat = resample(nfdat,fdat.shape[3],axis=3)
    nnfdat = resample(nnfdat,fdat.shape[3],axis=3)
    
    ofmridat = nib.load(os.path.join(datpath,func)) 
    
    funcdat = nib.Nifti1Image(nfdat,ofmridat.affine,dtype=ofmridat.get_data_dtype())
    nfuncdat = nib.Nifti1Image(nnfdat,ofmridat.affine,dtype=ofmridat.get_data_dtype())
    
    if isasl: 
        pref = 'cd'
        nboldendf = '_dune-asl_asl.nii'
    else: 
        pref = 'cd'
        nboldendf = '_dune-nonbold_bold.nii'

    if '_echo-' in func:
        splitf = func.split('_echo-')
        tfunc = splitf[0]+'_dune-bold_bold.nii'
        tnfunc = splitf[0]+nboldendf
    else:
        splitf = func.split('_bold')
        tfunc = splitf[0]+'_dune-bold_bold.nii'
        tnfunc = splitf[0]+nboldendf
        
    edfunc = os.path.join(datpath,pref+tfunc)
    nib.save(funcdat,edfunc)
    
    ednfunc = os.path.join(datpath,pref+tnfunc)
    nib.save(nfuncdat,ednfunc)
    
    return edfunc