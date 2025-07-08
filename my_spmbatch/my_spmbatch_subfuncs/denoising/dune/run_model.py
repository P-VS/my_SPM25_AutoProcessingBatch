#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 15:21:01 2024

@author: dr. Peter Van Schuerbeek
"""

import os
import utils
import torch_model

import torch
import torch.nn as nn
import torch.optim as optim 

import math
import numpy as np
import nibabel as nib

from scipy.stats.mstats import zscore
from copy import deepcopy

import json
 


def filt_signal(inSignal, Fs, Ny, startF,stopF):

    fft_inSignal = np.fft.rfft(inSignal,axis=2)
    
    # Determine which frequencies are included (assuming the rows range from 0Hz to Nyquist)
    f = Ny * (np.array(list(range(1, fft_inSignal.shape[2] + 1)))) / (fft_inSignal.shape[2])
    
    # Only include frequencies higher than 0.01Hz and lower than 0.12 Hz
    fexcl = np.squeeze(np.array(np.where(np.logical_or(f < startF,f > stopF))))
    
    fft_inSignal[:,:,fexcl] = 0
    
    filt_inSignal = np.fft.irfft(fft_inSignal,inSignal.shape[2],axis=2)
    
    return filt_inSignal
    

#-------------------------------------------------------------------------
 
def train_model(fdat,gmwmdata,nbdata,noise_file,jsonfile,datpath,load_trained_file,save_trained_file,isasl):
     
    train_gm, train_nb, train_noise = utils.get_train_data(fdat,gmwmdata,nbdata,noise_file,jsonfile)
    
    nvox = train_gm.shape[0]
    edim = train_gm.shape[1]
    tdim = train_gm.shape[2]

    TE=[]
    for ie in range(edim):
        if edim>1:
            splitf = jsonfile.split('_echo-1')
            tjsonfile = splitf[0]+'_echo-'+str(ie+1)+splitf[1]
        else:
            tjsonfile = jsonfile
            
        with open(tjsonfile,'r') as f: jsondat=json.load(f)
        f.close()
                        
        t_r = jsondat["RepetitionTime"]
        TE.append(jsondat["EchoTime"]*1000)

    # Determine sample frequency
    Fs = 1/t_r
    
    # Determine Nyquist-frequency
    Ny = Fs/2
    
    if isasl:
        train_gm = filt_signal(train_gm, Fs, Ny, 0.008, Ny-0.008)
        train_nb = filt_signal(train_nb, Fs, Ny, 0.008, Ny-0.008)
    
    """Train model"""
    t_train_gm = torch.from_numpy(train_gm)
    t_train_nb = torch.from_numpy(train_nb)
    t_train_noise = torch.from_numpy(train_noise)
    
    d_model = torch_model.DUNE(inChannels=edim,outlength=tdim)
    d_model.eval()
    
    try:
        d_model.load_state_dict(torch.load(load_trained_file, weights_only=False))
        num_epochs=30
    except:
        print('No valid trained_dune_model.pth file found.')
        num_epochs=30
        
    d_model.eval
    
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu") 
    d_model.to(device) 
    
    loss_fn = torch_model.custom_loss(TE,Ny)
    optimizer = optim.Adam(d_model.parameters(),lr=0.001,betas=(0.001,0.999),eps=1e-08,weight_decay=0)
    
    batch_size=5000
    
    # Split the training set into training and validation sets 
    val_percent = 0.2
    val_size    = int(val_percent * nvox)
    train_size  = nvox - val_size 
    
    train_gm_dataset, val_gm_dataset = torch.utils.data.random_split(t_train_gm,[train_size,val_size]) 
    train_nb_dataset, val_nb_dataset = torch.utils.data.random_split(t_train_nb,[train_size,val_size]) 
    train_noise_dataset, val_noise_dataset = torch.utils.data.random_split(t_train_noise,[train_size,val_size]) 
    
    # Train the model 
    prev_loss1 = 100.0
    prev_loss2 = 100.0
    prev_loss3 = 100.0
    
    prev_val_loss1 = 100.0
    prev_val_loss2 = 100.0
    prev_val_loss3 = 100.0
    save_loss = 100.0
 
    for epoch in range(num_epochs):  # loop over the dataset multiple times  
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
                    
            # zero the parameter gradients
            optimizer.zero_grad()
            
            outputs_gm,fft_outputs_gm = d_model(inputs_gm)
            
            with torch.no_grad():
                outputs_nb,_ = d_model(inputs_nb)
                
            loss = loss_fn(outputs_gm,outputs_nb,inputs_noise,inputs_gm,inputs_nb,fft_outputs_gm)
            
            optimizer.zero_grad() 
            loss.backward()
            optimizer.step()
            
            running_loss = loss.item()
        
            print('\rEpoch [{}/{}], Batch [{}/{}]: loss... {:.4f}'.format(epoch+1,num_epochs,i+1,len(train_gm_loader),running_loss), end="")
        
        # Iterate over data.
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
                outputs_nb,_ = d_model(inputs_nb)
                
                loss = loss_fn(outputs_gm,outputs_nb,inputs_noise,inputs_gm,inputs_nb,fft_outputs_gm)
                optimizer.step()
                
                running_val_loss = loss.item()
                
        print('\rEpoch [{}/{}]: Loss={:.4f}, Validation Loss={:.4f}\n'.format(epoch+1,num_epochs,running_loss,running_val_loss))
        
        if (running_loss+running_val_loss)<save_loss:
            torch.save(d_model.state_dict(),save_trained_file)
            
            save_loss = running_loss+running_val_loss

        if epoch>2:
            if (prev_loss1>running_loss) and (prev_val_loss1>running_val_loss):
                div1_loss= prev_loss1 - running_loss
                div1_vall_loss= prev_val_loss1 - running_val_loss
                
                if (div1_loss<0.01) and (div1_vall_loss<0.01): break
            
            if (running_val_loss>prev_val_loss1) and (prev_val_loss1>prev_val_loss2) and (prev_val_loss2>prev_val_loss3): break
            
            max_loss = max(running_loss,prev_loss1,prev_loss2,prev_loss3)
            min_loss = min(running_loss,prev_loss1,prev_loss2,prev_loss3)
            
            if (max_loss - min_loss)<0.02: break
        
        prev_loss3 = prev_loss2
        prev_loss2 = prev_loss1
        prev_loss1 = running_loss
        
        prev_val_loss3 = prev_val_loss2
        prev_val_loss2 = prev_val_loss1
        prev_val_loss1 = running_val_loss
        
    d_model.load_state_dict(torch.load(save_trained_file, weights_only=False))
    d_model.eval
            
    return d_model
    
#-------------------------------------------------------------------------
 
def apply_model(fdat,maskdata,d_model,func,datpath,jsonfile,isasl):
    
    func_data = utils.get_apply_func_data(fdat,maskdata)

    edim = func_data.shape[1]
    
    with open(jsonfile,'r') as f: jsondat=json.load(f)
    f.close()
                   
    t_r = jsondat["RepetitionTime"]

    # Determine sample frequency
    Fs = 1/t_r
   
    # Determine Nyquist-frequency
    Ny = Fs/2
    
    if isasl: 
        label_data = filt_signal(func_data, Fs, Ny, Ny-0.008, Ny)
        func_data = filt_signal(func_data, Fs, Ny, 0.008, Ny-0.008)
         
    dbold_data = np.zeros((func_data.shape[0],func_data.shape[2]))
    dnbold_data = np.zeros((func_data.shape[0],func_data.shape[2]))
    
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu") 
    d_model.to(device) 
    
    nsteps = math.ceil(func_data.shape[0]/5000)
    with torch.no_grad():
        for i in range(nsteps):
            last_vox = min((i+1)*5000,func_data.shape[0])
            
            fbatch = func_data[i*5000:last_vox,:,:]
            
            t_fbatch = torch.from_numpy(fbatch)
            inputs_gm = t_fbatch.float().to(device)
            
            m_output_fmri,_=d_model(inputs_gm)
            output_fmri = m_output_fmri[:,0,:]
            output_nfmri = m_output_fmri[:,1,:]
            
            dbold_data[i*5000:last_vox,:] = np.squeeze(output_fmri.numpy())
            dnbold_data[i*5000:last_vox,:] = np.squeeze(output_nfmri.numpy())
            
            print('\rBatch [{},{}]'.format(i+1,nsteps), end="")
            
    if isasl: dnbold_data = dnbold_data +  label_data[:,0,:]
    
    print('\n Saving model result \n')
    
    dfunf_file = utils.save_denoised_data(dbold_data,dnbold_data,fdat,maskdata,func,datpath,isasl)
    
    return dfunf_file
    
#-------------------------------------------------------------------------
 
def load_trained_model(edim,fdat,trained_file):
    
    tdim = fdat.shape[0]
    
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu") 
    
    d_model = torch_model.DUNE(inChannels=edim,outlength=tdim)
    d_model.eval()
    
    d_model.load_state_dict(torch.load(trained_file, weights_only=False))
    d_model.eval
    d_model.to(device) 
    
    return d_model