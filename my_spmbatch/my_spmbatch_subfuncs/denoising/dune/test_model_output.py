#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 15:23:59 2024

@author: dr. Peter Van Schuerbeek
"""

import os
import utils
import torch_model
import run_model
import json
import random

import torch
import torch.nn as nn

import math
import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt
from scipy.stats.mstats import zscore

def dunnet_model_performance(d_model,train_gm,train_nb,train_noise,jsonfile,isasl):
    
    nvox = train_gm.shape[0]
    edim = train_gm.shape[1]
    tdim = train_gm.shape[2]
    
    voxel = random.randint(0,nvox)
           
    with open(jsonfile,'r') as f: jsondat=json.load(f)
    f.close()
                    
    t_r = jsondat["RepetitionTime"]
    
    # Determine sample frequency
    Fs = 1/t_r
    
    # Determine Nyquist-frequency
    Ny = Fs/2
    
    if isasl:
        stopF = Ny-0.008
    else:
        stopF = Ny
    
    train_gm = run_model.filt_signal(train_gm, Fs, Ny, 0.008, stopF)
    train_nb = run_model.filt_signal(train_nb, Fs, Ny, 0.008, stopF)
       
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    d_model.to(device)
    
    ftrain_gm = train_gm.copy()
    
    t_train_gm = torch.from_numpy(train_gm)
    inputs_gm = t_train_gm.float().to(device)
    
    with torch.no_grad():
        output_fmri,_ = d_model(inputs_gm)
        
    denoised_fmri=output_fmri.numpy()
    
    train_gm = train_gm - np.mean(train_gm,axis=2,keepdims=True)
    train_gm = np.nan_to_num(train_gm / np.std(train_gm,axis=2,keepdims=True),nan=0.0)
    
    #Correlation denoised brain signal and denoised non-brain signals
    t_train_nb = torch.from_numpy(train_nb)
    inputs_nb = t_train_nb.float().to(device)
    
    with torch.no_grad():
        output_nb,_=d_model(inputs_nb)
        
    denoised_nb=output_nb.numpy()
    
    denoised_nb = denoised_nb - np.mean(denoised_nb,axis=2,keepdims=True)
    denoised_nb = np.nan_to_num(denoised_nb / np.std(denoised_nb,axis=2,keepdims=True),nan=0.0)
    
    train_nb = train_nb - np.mean(train_nb,axis=2,keepdims=True)
    train_nb = np.nan_to_num(train_nb / np.std(train_nb,axis=2,keepdims=True),nan=0.0)
    
    plt.figure()
    plt.plot(list(range(train_gm.shape[2])),denoised_fmri[voxel,1,:],color='g',linestyle='--',linewidth=0.8)
    plt.plot(list(range(train_gm.shape[2])),denoised_fmri[voxel,0,:],color='b',linestyle='--',linewidth=1.0)
    plt.title('Denoised brain signal')
    plt.legend(['Denoised non-BOLD signal','Denoised BOLD signal'],loc="upper right") 

    denoised_fmri = denoised_fmri - np.mean(denoised_fmri,axis=2,keepdims=True)
    denoised_fmri = np.nan_to_num(denoised_fmri / np.std(denoised_fmri,axis=2,keepdims=True),nan=0.0)
    
    plt.figure()
    plt.plot(list(range(train_gm.shape[2])),train_gm[voxel,0,:],color='y',linestyle='--',linewidth=0.6)
    plt.plot(list(range(train_gm.shape[2])),train_gm[voxel,1,:],color='m',linestyle='--',linewidth=0.8)
    if edim>2:
        plt.plot(list(range(train_gm.shape[2])),train_gm[voxel,2,:],color='g',linestyle='--',linewidth=0.6)
        plt.plot(list(range(train_gm.shape[2])),train_gm[voxel,3,:],color='c',linestyle='--',linewidth=0.6)
        plt.legend(['Echo 1','Echo 2','Echo 3','Echo 4','Denoised BOLD signal'],loc="lower right") 
    else:    
        plt.legend(['Echo 1','Echo 2','Denoised BOLD signal'],loc="lower right") 
    plt.plot(list(range(train_gm.shape[2])),denoised_fmri[voxel,0,:],color='b',linestyle='--',linewidth=1.2)
    plt.title('Denoised versus original signals')
       
    plt.figure()
    mtrain_gm = np.mean(train_gm,axis=1,keepdims=True)
    mtrain_gm = mtrain_gm - np.mean(mtrain_gm,axis=2,keepdims=True)
    mtrain_gm = np.nan_to_num(mtrain_gm / np.std(mtrain_gm,axis=2,keepdims=True),nan=0.0)
    
    plt.plot(list(range(train_gm.shape[2])),mtrain_gm[voxel,0,:],color='r',linestyle='--',linewidth=0.8)
    plt.plot(list(range(train_gm.shape[2])),denoised_fmri[voxel,0,:],color='b',linestyle='--',linewidth=1.0)
    plt.title('Denoised versus averaged signal')
    plt.legend(['Averaged signal','Denoised BOLD signal'],loc="lower right") 
    
    plt.figure()
    plt.plot(list(range(train_gm.shape[2])),train_nb[voxel,0,:],color='r',linestyle='--',linewidth=0.8)
    plt.plot(list(range(train_gm.shape[2])),denoised_nb[voxel,1,:],color='b',linestyle='--',linewidth=1.0)
    plt.title('Denoised signa versus TE1')
    plt.legend(['Echo 1','Denoised non-BOLD signal'],loc="lower right") 
    
    #Correlation denoised brain signal and noise regressors
    
    ocor_fmri_noise = np.sum(train_gm[:,0,:] * train_noise[:,0,:],axis=1)/(train_gm.shape[2]-1)
    mcor_fmri_noise = np.sum(mtrain_gm[:,0,:] * train_noise[:,0,:],axis=1)/(train_gm.shape[2]-1)
    cor_fmri_noise = np.sum(denoised_fmri[:,0,:] * train_noise[:,0,:],axis=1)/(train_gm.shape[2]-1)
    
    print('Mean correlation non-denoised fMRI and noise = {:.4f} ({:.4f})'.format(np.mean(ocor_fmri_noise),np.std(ocor_fmri_noise)))
    print('Mean correlation averaged fMRI and noise = {:.4f} ({:.4f})'.format(np.mean(mcor_fmri_noise),np.std(mcor_fmri_noise)))
    print('Mean correlation denoised fMRI and noise = {:.4f} ({:.4f})'.format(np.mean(cor_fmri_noise),np.std(cor_fmri_noise)))
    print('')
    
    plt.figure()
    plt.hist(ocor_fmri_noise, bins=200, range=([-1,1]), facecolor='r', alpha=0.5)
    plt.hist(cor_fmri_noise, bins=200, range=([-1,1]), facecolor='b', alpha=0.5)
    plt.xlim([-1.0,1.0])
    plt.title('Correlations between denoised fMRI and noise')
    plt.legend(['Without denoising','Denoised'],loc="upper right") 
    plt.axvline(x = 0, color = 'k')
    
    plt.figure()
    plt.hist(ocor_fmri_noise, bins=200, range=([-1,1]), facecolor='r', alpha=0.5)
    plt.hist(mcor_fmri_noise, bins=200, range=([-1,1]), facecolor='g', alpha=0.5)
    plt.hist(cor_fmri_noise, bins=200, range=([-1,1]), facecolor='b', alpha=0.5)
    plt.xlim([-1.0,1.0])
    plt.title('Correlations between brain signals and noise')
    plt.legend(['Without denoising','Averaged signal','Denoised'],loc="upper right") 
    plt.axvline(x = 0, color = 'k')

    mtrain_nb = np.mean(train_nb,axis=1,keepdims=True)
    
    ocor_fmri_nb = np.sum(train_gm[:,0,:] * train_nb[:,0,:],axis=1)/(train_gm.shape[2]-1)
    mcor_fmri_nb = np.sum(mtrain_gm[:,0,:] * mtrain_nb[:,0,:],axis=1)/(train_gm.shape[2]-1)
    cor_fmri_nb = np.sum(denoised_fmri[:,0,:] * denoised_nb[:,0,:],axis=1)/(train_gm.shape[2]-1)
    
    print('Mean correlation non-denoised fMRI and non-brain signals = {:.4f} ({:.4f})'.format(np.mean(ocor_fmri_nb),np.std(ocor_fmri_nb)))
    print('Mean correlation averaged brain and non-brain signals = {:.4f} ({:.4f})'.format(np.mean(mcor_fmri_nb),np.std(mcor_fmri_nb)))
    print('Mean correlation denoised BOLD in brain and non-brain voxels = {:.4f} ({:.4f})'.format(np.mean(cor_fmri_nb),np.std(cor_fmri_nb)))
    print('')
    
    plt.figure()
    plt.hist(ocor_fmri_nb, bins=200, range=([-1,1]), facecolor='r', alpha=0.5)
    plt.hist(cor_fmri_nb, bins=200, range=([-1,1]), facecolor='b', alpha=0.5)
    plt.xlim([-1.0,1.0])
    plt.title('Correlations between brain fMRI and non-brain signals')
    plt.legend(['Without denoising','Denoised'],loc="upper right") 
    plt.axvline(x = 0, color = 'k')
    
    plt.figure()
    plt.hist(ocor_fmri_nb, bins=200, range=([-1,1]), facecolor='r', alpha=0.5)
    plt.hist(mcor_fmri_nb, bins=200, range=([-1,1]), facecolor='g', alpha=0.5)
    plt.hist(cor_fmri_nb, bins=200, range=([-1,1]), facecolor='b', alpha=0.5)
    plt.xlim([-1.0,1.0])
    plt.title('Correlations between brain fMRI and non-brain signals')
    plt.legend(['Without denoising','Averaged signal','Denoised'],loc="upper right") 
    plt.axvline(x = 0, color = 'k')
    
    #Correlation between denoised brain signal
    
    train_gm2 = train_gm.copy()
    idx = torch.randperm(t_train_gm.shape[2])
    train_gm2 = train_gm2[:,:,idx]
    
    denoised_fmri2=output_fmri.numpy()
    idx = torch.randperm(output_fmri.shape[2])
    denoised_fmri2 = denoised_fmri2[:,:,idx]
    
    mtrain_gm2=mtrain_gm.copy()
    idx = torch.randperm(mtrain_gm.shape[2])
    mtrain_gm2 = mtrain_gm2[:,:,idx]
    
    ocor_fmri12 = np.sum(train_gm[:,0,:] * train_gm2[:,0,:],axis=1)/(train_gm.shape[2]-1)
    mcor_fmri12 = np.sum(mtrain_gm[:,0,:] * mtrain_gm2[:,0,:],axis=1)/(train_gm.shape[2]-1)
    cor_fmri12 = np.sum(denoised_fmri[:,0,:] * denoised_fmri2[:,0,:],axis=1)/(train_gm.shape[2]-1)
    
    print('Mean correlation between non-denoised brain signals = {:.4f} ({:.4f})'.format(np.mean(ocor_fmri12),np.std(ocor_fmri12)))
    print('Mean correlation between averaged brain signals = {:.4f} ({:.4f})'.format(np.mean(mcor_fmri12),np.std(mcor_fmri12)))
    print('Mean correlation denoised brain signals = {:.4f} ({:.4f})'.format(np.mean(cor_fmri12),np.std(cor_fmri12)))
    print('')
    
    plt.figure()
    plt.hist(ocor_fmri12, bins=200, range=([-1,1]), facecolor='r', alpha=0.5)
    plt.hist(cor_fmri12, bins=200, range=([-1,1]), facecolor='b', alpha=0.5)
    plt.xlim([-1.0,1.0])
    plt.title('Correlations between brain signals')
    plt.legend(['Without denoising','Denoised'],loc="upper right") 
    plt.axvline(x = 0, color = 'k')
    
    plt.figure()
    plt.hist(ocor_fmri12, bins=200, range=([-1,1]), facecolor='r', alpha=0.5)
    plt.hist(mcor_fmri12, bins=200, range=([-1,1]), facecolor='g', alpha=0.5)
    plt.hist(cor_fmri12, bins=200, range=([-1,1]), facecolor='b', alpha=0.5)
    plt.xlim([-1.0,1.0])
    plt.title('Correlations between brain signals')
    plt.legend(['Without denoising','Averaged signalls','Denoised'],loc="upper right") 
    plt.axvline(x = 0, color = 'k')
    
    return 1
    
#-------------------------------------------------------------------------
 
def loss_test(d_model,train_gm,train_nb,train_noise,jsonfile):
    
    nvox = train_gm.shape[0]
    edim = train_gm.shape[1]
    tdim = train_gm.shape[2]
    
    t_train_gm = torch.from_numpy(train_gm)
    t_train_nb = torch.from_numpy(train_nb)
    t_train_noise = torch.from_numpy(train_noise)
    
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    d_model.to(device)
    
    with open(jsonfile,'r') as f: jsondat=json.load(f)
    f.close()
                    
    t_r = jsondat["RepetitionTime"]
    
    # Determine sample frequency
    Fs = 1/t_r
    
    # Determine Nyquist-frequency
    Ny = Fs/2
    
    loss_fn = torch_model.custom_loss(Ny)
    
    # Split the training set into training and validation sets 
    batch_size=5000
    
    val_percent = 0.2
    val_size    = int(val_percent * nvox)
    train_size  = nvox - val_size 
    
    train_gm_dataset, val_gm_dataset = torch.utils.data.random_split(t_train_gm,[train_size,val_size]) 
    train_nb_dataset, val_nb_dataset = torch.utils.data.random_split(t_train_nb,[train_size,val_size]) 
    train_noise_dataset, val_noise_dataset = torch.utils.data.random_split(t_train_noise,[train_size,val_size]) 
    
    train_gm_loader = torch.utils.data.DataLoader(train_gm_dataset,batch_size=batch_size,shuffle=True,pin_memory=True) 
    val_gm_loader = torch.utils.data.DataLoader(val_gm_dataset,batch_size=batch_size,shuffle=False,pin_memory=True) 
    
    train_nb_loader = torch.utils.data.DataLoader(train_nb_dataset,batch_size=batch_size,shuffle=True,pin_memory=True) 
    val_nb_loader = torch.utils.data.DataLoader(val_nb_dataset,batch_size=batch_size,shuffle=False,pin_memory=True) 
    
    train_noise_loader = torch.utils.data.DataLoader(train_noise_dataset,batch_size=batch_size,shuffle=True,pin_memory=True) 
    val_noise_loader = torch.utils.data.DataLoader(val_noise_dataset,batch_size=batch_size,shuffle=False,pin_memory=True) 
  
    # Iterate over data.
    dataloader_gm_iter_in = iter(train_gm_loader)
    dataloader_nb_iter_in = iter(train_nb_loader)
    dataloader_noise_iter_in = iter(train_noise_loader)
    for i in range(len(train_gm_loader)):
        inputs_gm = next(dataloader_gm_iter_in)
        inputs_nb = next(dataloader_nb_iter_in)
        inputs_noise = next(dataloader_noise_iter_in)
        
        inputs_gm = inputs_gm.float().to(device)
        inputs_nb = inputs_nb.float().to(device)
        inputs_noise = inputs_noise.float().to(device)
        
        with torch.no_grad():
            outputs_gm, fft_outputs_gm = d_model(inputs_gm)
            outputs_nb, fft_outputs_nb = d_model(inputs_nb)
            
        loss = loss_fn(outputs_gm,outputs_nb,inputs_noise,inputs_gm,inputs_nb,fft_outputs_gm)
        
        running_loss = loss.item()
        
    dataloader_gm_iter_in = iter(val_gm_loader)
    dataloader_nb_iter_in = iter(val_nb_loader)
    dataloader_noise_iter_in = iter(val_noise_loader)
    with torch.no_grad():
        for i in range(len(val_gm_loader)):
            inputs_gm = next(dataloader_gm_iter_in)
            inputs_nb = next(dataloader_nb_iter_in)
            inputs_noise = next(dataloader_noise_iter_in)
            
            inputs_gm = inputs_gm.float().to(device)
            inputs_nb = inputs_nb.float().to(device)
            inputs_noise = inputs_noise.float().to(device)
            
            outputs_gm,fft_outputs_gm = d_model(inputs_gm)
            outputs_nb,fft_outputs_nb = d_model(inputs_nb)
            
            loss = loss_fn(outputs_gm,outputs_nb,inputs_noise,inputs_gm,inputs_nb,fft_outputs_gm)
            
            running_val_loss = loss.item()
            
    print('Loss={:.4f}, Validation Loss={:.4f}\n'.format(running_loss,running_val_loss))
    
    return 1
    
#-------------------------------------------------------------------------
 
class layers_DUNE(nn.Module):
    def __init__(self, in_channels,outlength):
        super(layers_DUNE,self).__init__()
        self.in_channels = in_channels
        self.out_channels = 1
        
        self.indouble = torch_model.doubleConv(in_channels, 16, outlength, kernel_size=7)
        self.down1 = torch_model.Encoder(16, 32, int(outlength/2)+1, kernel_size=5)
        self.down2 = torch_model.Encoder(32, 64, int(outlength/4)+1, kernel_size=3)
        
        self.middledouble = torch_model.doubleConv(64, 64, int(outlength/4)+1, kernel_size=3)
        
        self.up1 = torch_model.Decoder(64, 32, int(outlength/4)+1, kernel_size=3)        
        self.up2 = torch_model.Decoder(64, 16, int(outlength/2)+1, kernel_size=3)
        self.outdouble = torch_model.doubleConv(32, 2, outlength, kernel_size=1)
        
    def forward(self, x):
        
        x1 = self.indouble(x)
        x2 = self.down1(x1)
        x3 = self.down2(x2)
        
        x4 = self.middledouble(x3)
             
        x5 = self.up1(x4)
        x6 = self.up2(torch.cat([x5,x2],1))
        x7 = self.outdouble(torch.cat([x6,x1],1)) 
        
        #xf = fft.rfft(x,dim=2)
        
        return x1, x2, x3, x4, x5, x6, x7

#-------------------------------------------------------------------------
 
def layers_output_test(train_gm,datpath):
    
    nvox = train_gm.shape[0]
    edim = train_gm.shape[1]
    tdim = train_gm.shape[2]
    
    t_train_gm = torch.from_numpy(train_gm[0:10,:,:])
    
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    d_model = layers_DUNE(in_channels=edim,outlength=tdim)
    d_model.eval()
    
    dune_path = os.path.dirname(__file__)
    dune_file = os.path.join(dune_path,'trained_dune_model.pth')
    
    d_model.load_state_dict(torch.load(dune_file, weights_only=False))
    d_model.eval
    d_model.to(device) 
    
    inputs_gm = t_train_gm.float().to(device)
    with torch.no_grad():
        x1, x2, x3, x4, x5, x6, x7 = d_model(inputs_gm)
        
    x1=x1.numpy()
    x2=x2.numpy()
    x3=x3.numpy()
    x4=x4.numpy()
    x5=x5.numpy()
    x6=x6.numpy()
    x7=x7.numpy()
    
    plt.figure()
    plt.plot(list(range(train_gm.shape[2])),train_gm[1,0,:],color='r',linestyle='--',linewidth=0.8)
    plt.plot(list(range(train_gm.shape[2])),train_gm[1,1,:],color='b',linestyle='--',linewidth=0.8)    
    plt.title('original signal')
    plt.legend(['TE(1)','TE(2)'],loc="lower right") 
    plt.savefig(os.path.join(datpath,'OriginalSignal.png'))
    plt.close()
    
    ngraphs = 16
    ncols = 4
    nrows = 4
                    
    nplots = 1
    for ip in range(nplots):
        fig, axis = plt.subplots(nrows,ncols)
        fig.suptitle('Output of DUNE initial level', fontsize=10)
        
        for ir in range(nrows):
            for ic in range(ncols):
                lout = x1[1,ip*ngraphs+ir*ncols+ic,:]
                axis[ir, ic].plot(range(lout.shape[0]), lout,color='k',linestyle='--',linewidth=0.8) 
                axis[ir, ic].tick_params(axis='both', labelsize=5)
                
        fig.savefig(os.path.join(datpath,'OutputInitialLevel_'+str(ip)+'.png'))
        plt.close(fig=fig)
        
    nplots = 2
    for ip in range(nplots):
        fig, axis = plt.subplots(nrows,ncols)
        fig.suptitle('Output of DUNE down level 1', fontsize=10)
        
        for ir in range(nrows):
            for ic in range(ncols):
                lout = x2[1,ip*ngraphs+ir*ncols+ic,:]
                axis[ir, ic].plot(range(lout.shape[0]), lout,color='k',linestyle='--',linewidth=0.8) 
                axis[ir, ic].tick_params(axis='both', labelsize=5)
                
        fig.savefig(os.path.join(datpath,'OutputDownLevel1_'+str(ip)+'.png'))
        plt.close(fig=fig)
         
    nplots = 4
    for ip in range(nplots):
        fig, axis = plt.subplots(nrows,ncols)
        fig.suptitle('Output of DUNE down level 2', fontsize=10)
        
        for ir in range(nrows):
            for ic in range(ncols):
                lout = x3[1,ip*ngraphs+ir*ncols+ic,:]
                axis[ir, ic].plot(range(lout.shape[0]), lout,color='k',linestyle='--',linewidth=0.8) 
                axis[ir, ic].tick_params(axis='both', labelsize=5)
                
        fig.savefig(os.path.join(datpath,'OutputDownLevel2_'+str(ip)+'.png'))
        plt.close(fig=fig)
        
    nplots = 4
    for ip in range(nplots):
        fig, axis = plt.subplots(nrows,ncols)
        fig.suptitle('Output of DUNE middle level 2', fontsize=10)
        
        for ir in range(nrows):
            for ic in range(ncols):
                lout = x4[1,ip*ngraphs+ir*ncols+ic,:]
                axis[ir, ic].plot(range(lout.shape[0]), lout,color='k',linestyle='--',linewidth=0.8) 
                axis[ir, ic].tick_params(axis='both', labelsize=5)
                
        fig.savefig(os.path.join(datpath,'OutputMiddleLevel2_'+str(ip)+'.png'))
        plt.close(fig=fig)
        
    nplots = 2
    for ip in range(nplots):
        fig, axis = plt.subplots(nrows,ncols)
        fig.suptitle('Output of DUNE Up level 2', fontsize=10)
        
        for ir in range(nrows):
            for ic in range(ncols):
                lout = x5[1,ip*ngraphs+ir*ncols+ic,:]
                axis[ir, ic].plot(range(lout.shape[0]), lout,color='k',linestyle='--',linewidth=0.8) 
                axis[ir, ic].tick_params(axis='both', labelsize=5)
                
        fig.savefig(os.path.join(datpath,'OutputUpLevel2_'+str(ip)+'.png'))
        plt.close(fig=fig)
        
    nplots = 1
    for ip in range(nplots):
        fig, axis = plt.subplots(nrows,ncols)
        fig.suptitle('Output of DUNE Up level 1', fontsize=10)
        
        for ir in range(nrows):
            for ic in range(ncols):
                lout = x6[1,ip*ngraphs+ir*ncols+ic,:]
                axis[ir, ic].plot(range(lout.shape[0]), lout,color='k',linestyle='--',linewidth=0.8) 
                axis[ir, ic].tick_params(axis='both', labelsize=5)
                
        fig.savefig(os.path.join(datpath,'OutputUpLevel1_'+str(ip)+'.png'))
        plt.close(fig=fig)
        
    plt.figure()
    plt.plot(list(range(train_gm.shape[2])),x7[1,1,:],color='r',linestyle='--',linewidth=0.8)
    plt.plot(list(range(train_gm.shape[2])),x7[1,0,:],color='b',linestyle='--',linewidth=0.8)    
    plt.title('Output signal')
    plt.legend(['non-BOLD','TBOLD'],loc="lower right") 
    plt.savefig(os.path.join(datpath,'OutputDUNE.png'))
    plt.close()
    
    return 1