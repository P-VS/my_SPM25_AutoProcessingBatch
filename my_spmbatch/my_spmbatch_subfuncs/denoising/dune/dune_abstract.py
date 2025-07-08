#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 11:15:01 2024

@author: dr. Peter Van Schuerbeek
"""

import os
import utils
import run_model
import json

import random
import math
import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt

from scipy.stats.mstats import zscore
    
def get_brain_and_non_brain_masks(gmfile,wmfile,csffile,maskfile):
    
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
    nbdata[csf_data>0.2]=1
    nbdata[gm_data+wm_data>0.1]=0
    nbdata[nbdata>0.0]=1
    
    print('make masks done')
    
    return gmwmdata, nbdata

"""
---------------------------------------------------------------------------------------
""" 
def get_func_data(funcfile,gmwmdata,nbdata):
    
    """
    Load brain (GM+WM) and no-brain data 
    """
    
    print('Start loading '+funcfile)

    fmridat = nib.load(funcfile)     
    fdat = fmridat.get_fdata()
    
    fdat = np.transpose(fdat,axes=(3,0,1,2))

    #gmfdat = fdat[:,gmwmdata>0.0].T
    all_nbfdat = fdat[:,nbdata>0.0].T
    print('2')
    print(all_nbfdat.shape)
    nvoxel_eval = gmfdat.shape[0]
    #nvoxel_eval = min(nvoxel_gmfdat,50000)
    
    #eval_gmfdat = np.random.permutation(int(nvoxel_gmfdat))[:int(nvoxel_eval)]
    #gmfdat = all_gmfdat[eval_gmfdat,:,:]
    
    nvoxel_nbdat = all_nbfdat.shape[0]
    if nvoxel_nbdat<nvoxel_eval:
        all_nbfdat = np.repeat(all_nbfdat, math.ceil(nvoxel_eval/nvoxel_nbdat), axis=0)
        nvoxel_nbdat = all_nbfdat.shape[0]
        
    eval_nbfdat = np.random.permutation(int(nvoxel_nbdat))[:int(nvoxel_eval)]
    nbfdat = all_nbfdat[eval_nbfdat,:]
    
    print('get functional data done')
    
    return gmfdat, nbfdat

"""
---------------------------------------------------------------------------------------
""" 

#def dune_model_comparison(subject):
    
subject = 4

if subject<10:
    substring = '0'+str(subject)
else:
    substring = str(subject)

datpath = '/Volumes/LaCie/UZ_Brussel/ME_fMRI_GE/data/sub-'+substring+'/ses-001/'
dune_folder = datpath+'preproc_func_dunnet/'
meica_folder = datpath+'preproc_EmoFaces_func_meica/'
nodenoise_folder = datpath+'preproc_EmoFaces_func/'

funcfile = 'swacduresub-'+substring+'_task-ME-EmoFaces_bold.nii' 

dune_funcfile = dune_folder+funcfile
meica_funcfile = meica_folder+funcfile
nodenoise_funcfile = nodenoise_folder+'swacuresub-'+substring+'_task-ME-EmoFaces_bold.nii' 

jsonfile = datpath+'func/sub-'+substring+'_task-ME-EmoFaces_echo-1_bold.json' 
noise_file = dune_folder+'der_rp_esub-'+substring+'_task-ME-EmoFaces_bold.txt'

gmfile = datpath+'func/wc1uresub-'+substring+'_task-ME-EmoFaces_echo-1_bold.nii' 
wmfile = datpath+'func/wc2uresub-'+substring+'_task-ME-EmoFaces_echo-1_bold.nii' 
csffile = datpath+'func/wc3uresub-'+substring+'_task-ME-EmoFaces_echo-1_bold.nii' 

with open(jsonfile,'r') as f: jsondat=json.load(f)
f.close()
            
t_r = jsondat["RepetitionTime"]

maskfile = dune_folder+'mask_swacduresub-'+substring+'_task-ME-EmoFaces_bold.nii' 

#Get brain and non-brain mask
gmwmdata, nbdata = get_brain_and_non_brain_masks(gmfile,wmfile,csffile,maskfile)

#Get funnctional data
dune_gmfdat, dune_nbfdat = get_func_data(dune_funcfile,gmwmdata,nbdata)
meica_gmfdat, meica_nbfdat = get_func_data(meica_funcfile,gmwmdata,nbdata)
nodenoise_gmfdat, nodenoise_nbfdat = get_func_data(nodenoise_funcfile,gmwmdata,nbdata)

random_voxel = random.randint(1,dune_gmfdat.shape[0])

timeax = [k*t_r for k in range(0,nodenoise_gmfdat.shape[1])]

#Example signal
plt.figure()
plt.plot(timeax,zscore(nodenoise_gmfdat[random_voxel,:]),color='r',linestyle='--',linewidth=0.8)
plt.plot(timeax,zscore(meica_gmfdat[random_voxel,:]),color='g',linestyle='--',linewidth=0.8)
plt.plot(timeax,zscore(dune_gmfdat[random_voxel,:]),color='b',linestyle='--',linewidth=1.0)
plt.title('Signal from a random voxel in the brain')
plt.legend(['default SPM','ME-ICA-AROMA','ME-DUNE'],loc="lower right") 
plt.savefig(os.path.join(datpath,'signal_randomVox_default-meica-dune.png'))

plt.figure()
plt.plot(timeax,zscore(nodenoise_gmfdat[random_voxel,:]),color='r',linestyle='--',linewidth=0.8)
plt.plot(timeax,zscore(meica_gmfdat[random_voxel,:]),color='b',linestyle='--',linewidth=1.0)
plt.title('Signal from a random voxel in the brain')
plt.legend(['default SPM','ME-ICA-AROMA'],loc="lower right") 
plt.savefig(os.path.join(datpath,'signal_randomVox_default-meica.png'))

plt.figure()
plt.plot(timeax,zscore(nodenoise_gmfdat[random_voxel,:]),color='r',linestyle='--',linewidth=0.8)
plt.plot(timeax,zscore(dune_gmfdat[random_voxel,:]),color='b',linestyle='--',linewidth=1.0)
plt.title('Signal from a random voxel in the brain')
plt.legend(['default SPM','ME-DUNE'],loc="lower right") 
plt.savefig(os.path.join(datpath,'signal_randomVox_default-dune.png'))

#tSNR
dune_tsnr_fdat = np.mean(dune_gmfdat,axis=1)
meica_tsnr_fdat = np.mean(meica_gmfdat,axis=1)
nodenoise_tsnr_fdat = np.mean(nodenoise_gmfdat,axis=1)

dune_tsnr_fdat = dune_tsnr_fdat / np.std(dune_gmfdat,axis=1)
meica_tsnr_fdat = meica_tsnr_fdat / np.std(meica_gmfdat,axis=1)
nodenoise_tsnr_fdat = nodenoise_tsnr_fdat / np.std(nodenoise_gmfdat,axis=1)

dune_tsnr_fdat = np.nan_to_num(dune_tsnr_fdat,nan=0.0)
meica_tsnr_fdat = np.nan_to_num(meica_tsnr_fdat,nan=0.0)
nodenoise_tsnr_fdat = np.nan_to_num(nodenoise_tsnr_fdat,nan=0.0)

print('')
print('Default SPM: tSNR= {:.2f} ({:.2f})'.format(np.mean(nodenoise_tsnr_fdat),np.std(nodenoise_tsnr_fdat)))
print('ME-ICA-AROMA: tSNR= {:.2f} ({:.2f})'.format(np.mean(meica_tsnr_fdat),np.std(meica_tsnr_fdat)))
print('ME-DUNE: tSNR= {:.2f} ({:.2f})'.format(np.mean(dune_tsnr_fdat),np.std(dune_tsnr_fdat)))
print('')

#Correlation brain / non-brain signals
dune_gmfdat = zscore(dune_gmfdat,axis=1)
meica_gmfdat = zscore(meica_gmfdat,axis=1)
nodenoise_gmfdat = zscore(nodenoise_gmfdat,axis=1)

dune_nbfdat = zscore(dune_nbfdat,axis=1)
meica_nbfdat = zscore(meica_nbfdat,axis=1)
nodenoise_nbfdat = zscore(nodenoise_nbfdat,axis=1)

cor_dune_gmnnb = np.sum(dune_gmfdat * dune_nbfdat,axis=1)/(dune_gmfdat.shape[1]-1)
cor_meica_gmnnb = np.sum(meica_gmfdat * meica_nbfdat,axis=1)/(meica_gmfdat.shape[1]-1)
cor_nodenoise_gmnnb = np.sum(nodenoise_gmfdat * nodenoise_nbfdat,axis=1)/(nodenoise_gmfdat.shape[1]-1)

print('')
print('Default SPM: correlation between brain and non-brain signals = {:.2f} ({:.2f})'.format(np.mean(cor_nodenoise_gmnnb),np.std(cor_nodenoise_gmnnb)))
print('ME-ICA-AROMA: correlation between brain and non-brain signals = {:.2f} ({:.2f})'.format(np.mean(cor_meica_gmnnb),np.std(cor_meica_gmnnb)))
print('ME-DUNE: correlation between brain and non-brain signals = {:.2f} ({:.2f})'.format(np.mean(cor_dune_gmnnb),np.std(cor_dune_gmnnb)))
print('')

plt.figure()
plt.hist(cor_nodenoise_gmnnb, bins=200, range=([-1,1]), facecolor='r', alpha=0.5)
plt.hist(cor_meica_gmnnb, bins=200, range=([-1,1]), facecolor='g', alpha=0.5)
plt.hist(cor_dune_gmnnb, bins=200, range=([-1,1]), facecolor='b', alpha=0.5)
plt.xlim([-1.0,1.0])
plt.title('Correlations between brain and non-brain signals')
plt.legend(['default SPM','ME-ICA-AROMA','ME-DUNE'],loc="upper right") 
plt.axvline(x = 0, color = 'k')
plt.savefig(os.path.join(datpath,'correlation_bnb_default-meica-dune.png'))

plt.figure()
plt.hist(cor_nodenoise_gmnnb, bins=200, range=([-1,1]), facecolor='r', alpha=0.5)
plt.hist(cor_meica_gmnnb, bins=200, range=([-1,1]), facecolor='g', alpha=0.5)
plt.xlim([-1.0,1.0])
plt.title('Correlations between brain and non-brain signals')
plt.legend(['default SPM','ME-ICA-AROMA'],loc="upper right") 
plt.axvline(x = 0, color = 'k')
plt.savefig(os.path.join(datpath,'correlation_bnb_default-meica.png'))

plt.figure()
plt.hist(cor_nodenoise_gmnnb, bins=200, range=([-1,1]), facecolor='r', alpha=0.5)
plt.hist(cor_dune_gmnnb, bins=200, range=([-1,1]), facecolor='b', alpha=0.5)
plt.xlim([-1.0,1.0])
plt.title('Correlations between brain and non-brain signals')
plt.legend(['default SPM','ME-DUNE'],loc="upper right") 
plt.axvline(x = 0, color = 'k')
plt.savefig(os.path.join(datpath,'correlation_bnb_default-dune.png'))

#Correlation brain / noise regressors
params = np.genfromtxt(noise_file)
all_noisedat = params.transpose()

nvoxel_noisedat = all_noisedat.shape[0]
if nvoxel_noisedat<dune_gmfdat.shape[0]:
    all_noisedat = np.repeat(all_noisedat, math.ceil(dune_gmfdat.shape[0]/nvoxel_noisedat), axis=0)
    nvoxel_noisedat = all_noisedat.shape[0]
    
indx_noisedat = np.random.permutation(int(nvoxel_noisedat))[:int(dune_gmfdat.shape[0])]    
noisedat = all_noisedat[indx_noisedat,:]

noisedat = zscore(noisedat,axis=1)

cor_dune_gmnoise = np.sum(dune_gmfdat * noisedat,axis=1)/(dune_gmfdat.shape[1]-1)
cor_meica_gmnoise = np.sum(meica_gmfdat * noisedat,axis=1)/(meica_gmfdat.shape[1]-1)
cor_nodenoise_gmnoise = np.sum(nodenoise_gmfdat * noisedat,axis=1)/(nodenoise_gmfdat.shape[1]-1)

print('')
print('Default SPM: correlation between brain and noise signals = {:.2f} ({:.2f})'.format(np.mean(cor_nodenoise_gmnoise),np.std(cor_nodenoise_gmnoise)))
print('ME-ICA-AROMA: correlation between brain and noise signals = {:.2f} ({:.2f})'.format(np.mean(cor_meica_gmnoise),np.std(cor_meica_gmnoise)))
print('ME-DUNE: correlation between brain and noise signals = {:.2f} ({:.2f})'.format(np.mean(cor_dune_gmnoise),np.std(cor_dune_gmnoise)))
print('')

plt.figure()
plt.hist(cor_nodenoise_gmnoise, bins=200, range=([-1,1]), facecolor='r', alpha=0.5)
plt.hist(cor_meica_gmnoise, bins=200, range=([-1,1]), facecolor='g', alpha=0.5)
plt.hist(cor_dune_gmnoise, bins=200, range=([-1,1]), facecolor='b', alpha=0.5)
plt.xlim([-1.0,1.0])
plt.title('Correlations between brain and noise signals')
plt.legend(['default SPM','ME-ICA-AROMA','ME-DUNE'],loc="upper right") 
plt.axvline(x = 0, color = 'k')
plt.savefig(os.path.join(datpath,'correlation_bnoise_default-meica-dune.png'))

plt.figure()
plt.hist(cor_nodenoise_gmnoise, bins=200, range=([-1,1]), facecolor='r', alpha=0.5)
plt.hist(cor_meica_gmnoise, bins=200, range=([-1,1]), facecolor='g', alpha=0.5)
plt.xlim([-1.0,1.0])
plt.title('Correlations between brain and noise signals')
plt.legend(['default SPM','ME-ICA-AROMA'],loc="upper right") 
plt.axvline(x = 0, color = 'k')
plt.savefig(os.path.join(datpath,'correlation_bnoise_default-meica.png'))

plt.figure()
plt.hist(cor_nodenoise_gmnoise, bins=200, range=([-1,1]), facecolor='r', alpha=0.5)
plt.hist(cor_dune_gmnoise, bins=200, range=([-1,1]), facecolor='b', alpha=0.5)
plt.xlim([-1.0,1.0])
plt.title('Correlations between brain and noise signals')
plt.legend(['default SPM','ME-DUNE'],loc="upper right") 
plt.axvline(x = 0, color = 'k')
plt.savefig(os.path.join(datpath,'correlation_bnoise_default-dune.png'))