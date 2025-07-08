#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 16:26:26 2024

@author: dr. Peter Van Schuerbeek (UZ Brussel)
"""

import sys
import utils
import torch_model
import run_model
import test_model_output

sub=[1,2,3,4,5,6,7,8,9,10,11,12,13]
    
for i in sub:
    print('Start DUNE for subject '+str(i))
    
    if i<10:
        substring = '0'+str(i)
    else:
        substring = str(i)
            
    datpath = '/Volumes/LaCie/UZ_Brussel/ASLBOLD_OpenNeuro_FT/IndData/sub-'+substring+'/ses-001/func/'
    funcfile_e1 = 'lafresub-'+substring+'_task-bilateralfingertapping_echo-1_aslbold.nii'
    jsonfile = '/Volumes/LaCie/UZ_Brussel/ASLBOLD_OpenNeuro_FT/IndData/sub-'+substring+'/ses-001/func/sub-'+substring+'_task-bilateralfingertapping_echo-1_aslbold.json'
    nechoes = 4
    
    c1file = '/Volumes/LaCie/UZ_Brussel/ASLBOLD_OpenNeuro_FT/IndData/sub-'+substring+'/ses-001/func/c1afresub-'+substring+'_task-bilateralfingertapping_echo-1_aslbold.nii'
    c2file = '/Volumes/LaCie/UZ_Brussel/ASLBOLD_OpenNeuro_FT/IndData/sub-'+substring+'/ses-001/func/c2afresub-'+substring+'_task-bilateralfingertapping_echo-1_aslbold.nii'
    c3file = '/Volumes/LaCie/UZ_Brussel/ASLBOLD_OpenNeuro_FT/IndData/sub-'+substring+'/ses-001/func/c3afresub-'+substring+'_task-bilateralfingertapping_echo-1_aslbold.nii'
    fmaskfile = '/Volumes/LaCie/UZ_Brussel/ASLBOLD_OpenNeuro_FT/IndData/sub-'+substring+'/ses-001/func/fmask_afresub-'+substring+'_task-bilateralfingertapping_echo-1_aslbold.nii'
    noise_file = '/Volumes/LaCie/UZ_Brussel/ASLBOLD_OpenNeuro_FT/IndData/sub-'+substring+'/ses-001/func/der_rp_esub-'+substring+'_task-bilateralfingertapping_aslbold.txt'
    
    load_trained_file = '/Volumes/LaCie/UZ_Brussel/ASLBOLD_OpenNeuro_FT/IndData/sub-'+substring+'/ses-001/func/trained_dune_model_05062025_2.pth'
    save_trained_file = '/Volumes/LaCie/UZ_Brussel/ASLBOLD_OpenNeuro_FT/IndData/sub-'+substring+'/ses-001/func/trained_dune_model_05062025_2.pth'
    isasl = True
    
    gmwmdata, nbdata, braindata = utils.get_brain_and_non_brain_masks(datpath,c1file,c2file,c3file,fmaskfile)
    
    fdat = utils.get_func_data(datpath,funcfile_e1,jsonfile,nechoes)
    
    d_model = run_model.train_model(fdat,gmwmdata,nbdata,noise_file,jsonfile,datpath,load_trained_file,save_trained_file,isasl)
    #d_model = run_model.load_trained_model(nechoes,fdat,save_trained_file)
    
    out = run_model.apply_model(fdat,braindata,d_model,funcfile_e1,datpath,jsonfile,isasl)
    
    #train_gm, train_nb, train_noise = utils.get_train_data(fdat,gmwmdata,nbdata,noise_file,jsonfile)
    #out=test_model_output.dunnet_model_performance(d_model,train_gm,train_nb,train_noise,jsonfile,isasl)
    #out = test_model_output.loss_test(d_model,train_gm,train_nb,train_noise,jsonfile)
    #out = test_model_output.layers_output_test(train_gm,datpath)
      
    """
    datpath = '/Volumes/LaCie/UZ_Brussel/ME_fMRI_GE/data/sub-'+substring+'/ses-001/func/'
    funcfile_e1 = 'fauresub-'+substring+'_task-ME-EmoFaces_echo-1_bold.nii'
    jsonfile = '/Volumes/LaCie/UZ_Brussel/ME_fMRI_GE/data/sub-'+substring+'/ses-001/func/sub-'+substring+'_task-ME-EmoFaces_echo-1_bold.json'
    nechoes = 2
    
    c1file = '/Volumes/LaCie/UZ_Brussel/ME_fMRI_GE/data/sub-'+substring+'/ses-001/func/c1auresub-'+substring+'_task-ME-EmoFaces_echo-1_bold.nii'
    c2file = '/Volumes/LaCie/UZ_Brussel/ME_fMRI_GE/data/sub-'+substring+'/ses-001/func/c2auresub-'+substring+'_task-ME-EmoFaces_echo-1_bold.nii'
    c3file = '/Volumes/LaCie/UZ_Brussel/ME_fMRI_GE/data/sub-'+substring+'/ses-001/func/c3auresub-'+substring+'_task-ME-EmoFaces_echo-1_bold.nii'
    fmaskfile = '/Volumes/LaCie/UZ_Brussel/ME_fMRI_GE/data/sub-'+substring+'/ses-001/func/fmask_auresub-'+substring+'_task-ME-EmoFaces_echo-1_bold.nii'
    noise_file = '/Volumes/LaCie/UZ_Brussel/ME_fMRI_GE/data/sub-'+substring+'/ses-001/func/acc_der_rp_esub-'+substring+'_task-ME-EmoFaces_bold.txt'
    
    load_trained_file = '/Volumes/LaCie/UZ_Brussel/ME_fMRI_GE/data/sub-'+substring+'/ses-001/func/trained_dune_model_WCOR-NBRAIN_ME-EmoFaces_20250519.pth'
    save_trained_file = '/Volumes/LaCie/UZ_Brussel/ME_fMRI_GE/data/sub-'+substring+'/ses-001/func/trained_dune_model_WCOR-NBRAIN_ME-EmoFaces_20250519.pth'
    isasl = False
    
    gmwmdata, nbdata, braindata = utils.get_brain_and_non_brain_masks(datpath,c1file,c2file,c3file,fmaskfile)
    
    fdat = utils.get_func_data(datpath,funcfile_e1,jsonfile,nechoes)
    
    d_model = run_model.train_model(fdat,gmwmdata,nbdata,noise_file,jsonfile,datpath,load_trained_file,save_trained_file,isasl)
    #d_model = run_model.load_trained_model(nechoes,fdat,save_trained_file)
    
    out = run_model.apply_model(fdat,braindata,d_model,funcfile_e1,datpath,jsonfile,isasl)
    
    #train_gm, train_nb, train_noise = utils.get_train_data(fdat,gmwmdata,nbdata,noise_file,jsonfile)
    #out=test_model_output.dunnet_model_performance(d_model,train_gm,train_nb,train_noise,jsonfile,isasl)
    #out = test_model_output.loss_test(d_model,train_gm,train_nb,train_noise,jsonfile)
    #out = test_model_output.layers_output_test(train_gm,datpath)
      
    datpath = '/Volumes/LaCie/UZ_Brussel/ME_fMRI_GE/data/sub-'+substring+'/ses-001/func/'
    funcfile_e1 = 'fauresub-'+substring+'_task-ME-EFT_echo-1_bold.nii'
    jsonfile = '/Volumes/LaCie/UZ_Brussel/ME_fMRI_GE/data/sub-'+substring+'/ses-001/func/sub-'+substring+'_task-ME-EFT_echo-1_bold.json'
    nechoes = 2
    
    c1file = '/Volumes/LaCie/UZ_Brussel/ME_fMRI_GE/data/sub-'+substring+'/ses-001/func/c1auresub-'+substring+'_task-ME-EFT_echo-1_bold.nii'
    c2file = '/Volumes/LaCie/UZ_Brussel/ME_fMRI_GE/data/sub-'+substring+'/ses-001/func/c2auresub-'+substring+'_task-ME-EFT_echo-1_bold.nii'
    c3file = '/Volumes/LaCie/UZ_Brussel/ME_fMRI_GE/data/sub-'+substring+'/ses-001/func/c3auresub-'+substring+'_task-ME-EFT_echo-1_bold.nii'
    fmaskfile = '/Volumes/LaCie/UZ_Brussel/ME_fMRI_GE/data/sub-'+substring+'/ses-001/func/fmask_auresub-'+substring+'_task-ME-EFT_echo-1_bold.nii'
    noise_file = '/Volumes/LaCie/UZ_Brussel/ME_fMRI_GE/data/sub-'+substring+'/ses-001/func/acc_der_rp_esub-'+substring+'_task-ME-EFT_bold.txt'
    
    load_trained_file = '/Volumes/LaCie/UZ_Brussel/ME_fMRI_GE/data/sub-'+substring+'/ses-001/func/trained_dune_model_WCOR-NBRAIN_ME-EFT_20250519.pth'
    save_trained_file = '/Volumes/LaCie/UZ_Brussel/ME_fMRI_GE/data/sub-'+substring+'/ses-001/func/trained_dune_model_WCOR-NBRAIN_ME-EFT_20250519.pth'
    isasl = False
    
    gmwmdata, nbdata, braindata = utils.get_brain_and_non_brain_masks(datpath,c1file,c2file,c3file,fmaskfile)
    
    fdat = utils.get_func_data(datpath,funcfile_e1,jsonfile,nechoes)
    
    d_model = run_model.train_model(fdat,gmwmdata,nbdata,noise_file,jsonfile,datpath,load_trained_file,save_trained_file,isasl)
    #d_model = run_model.load_trained_model(nechoes,fdat,save_trained_file)
    
    out = run_model.apply_model(fdat,braindata,d_model,funcfile_e1,datpath,jsonfile,isasl)
    """
    print('Done DUNE for subject '+str(i))
    
print('done')