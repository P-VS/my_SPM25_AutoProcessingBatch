  #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 15:26:35 2024

@author: dr. Peter Van Schuerbeek
"""

import torch
import torch.nn as nn
import torch.fft as fft

import math
import numpy as np

class doubleConv(nn.Module):
    def __init__(self, inChannels, outChannels, length, kernel_size=7):
        super().__init__()
        padding = int((kernel_size - 1) / 2)
        
        self.dConv1 = nn.Sequential(
            nn.Conv1d(inChannels,outChannels,kernel_size=kernel_size,bias=True,padding=padding),
            nn.LayerNorm((outChannels,length),bias=True)
        )
        self.dConv2 = nn.Sequential(
            nn.Conv1d(outChannels,outChannels,kernel_size=kernel_size,bias=True,padding=padding),
            nn.LayerNorm((outChannels,length),bias=True)
        )
            
    def forward(self, x):
        x = self.dConv1(x)
        x = self.dConv2(x)
        
        return x
    
#-------------------------------------------------------------------------

class Encoder(nn.Module):    
    def __init__(self, inChannels, outChannels, length, kernel_size=7):
        super().__init__()
        padding = int((kernel_size - 1) / 2)
        
        self.maxpool_conv = nn.Sequential(
            nn.MaxPool1d(kernel_size=kernel_size,stride=2,padding=padding),
            doubleConv(inChannels, outChannels, length, kernel_size)
        )
        
    def forward(self, x):
            
        return self.maxpool_conv(x)
    
class Decoder(nn.Module):    
    def __init__(self, inChannels, outChannels, length, kernel_size=7):
        super().__init__()
        padding = int((kernel_size - 1) / 2)
        
        if (length % 2):
            output_padding=0
        else:
            output_padding=1
        
        self.convUpS = nn.Sequential(
            doubleConv(inChannels, outChannels, length, kernel_size),
            nn.ConvTranspose1d(outChannels,outChannels,stride=2,bias=True,kernel_size=kernel_size,output_padding=output_padding,padding=padding),
            nn.LayerNorm((outChannels,(length * 2)-(length % 2)),bias=True)
        )
        
    def forward(self, x):
        x = self.convUpS(x)
            
        return x
   
#-------------------------------------------------------------------------
 
class DUNE(nn.Module):
    def __init__(self, inChannels,outlength):
        super(DUNE,self).__init__()
        self.inChannels = inChannels
        self.outChannels = 1
        
        self.indouble = doubleConv(inChannels, 16, outlength, kernel_size=7)
        self.down1 = Encoder(16, 32, int(outlength/2)+(outlength % 2), kernel_size=5)
        self.down2 = Encoder(32, 64, int(outlength/4)+(outlength % 2), kernel_size=3)
        
        self.middledouble = doubleConv(64, 64, int(outlength/4)+(outlength % 2), kernel_size=3)
        
        self.up1 = Decoder(64, 32, int(outlength/4)+(outlength % 2), kernel_size=3)        
        self.up2 = Decoder(64, 16, int(outlength/2)+(outlength % 2), kernel_size=3)   
        self.outdouble = doubleConv(32, 2, outlength, kernel_size=1)
        
    def forward(self, x):
        
        x1 = self.indouble(x)
        x2 = self.down1(x1)
        x3 = self.down2(x2)
        
        x = self.middledouble(x3)
         
        x = self.up1(x)
        x = self.up2(torch.cat([x,x2],1))
        x = self.outdouble(torch.cat([x,x1],1)) 
        
        x2 = x-torch.mean(x,dim=2,keepdim=True)
        x2 = x/torch.std(x2,dim=2,keepdim=True)
        xf = fft.rfft(x2,dim=2)
        
        return x, xf
    
#-------------------------------------------------------------------------
 
class custom_loss(nn.Module):
    def __init__(self, TE,Ny):
        super(custom_loss, self).__init__()
        
        self.Ny = Ny
        self.TE = TE
        
    def forward(self, y_pred, y_nbrain, y_noise, y_orig, y_nborig,fft_y_pred):
        
        """Minimal correlation between BOLD part and noise regressors"""  
        
        output_fMRI = y_pred[:,0,:]
        input_noise  = y_noise[:,0,:]
        
        output_fMRI = output_fMRI - torch.mean(output_fMRI,dim=1,keepdim=True)
        output_fMRI = output_fMRI / torch.std(output_fMRI,dim=1,keepdim=True)
        
        corr_output = torch.sum(output_fMRI*input_noise,dim=1)/(output_fMRI.shape[1]-1)
        
        loss_corr_fMRInoise = abs(torch.mean(torch.abs(corr_output))) 
        
        """Minimal correlation between non-BOLD part and noise regressors"""  
        
        output_S0 = y_pred[:,1,:]
        input_noise  = y_noise[:,0,:]
        
        output_S0 = output_S0 - torch.mean(output_S0,dim=1,keepdim=True)
        output_S0 = output_S0 / torch.std(output_S0,dim=1,keepdim=True)
        
        corr_output = torch.sum(output_S0*input_noise,dim=1)/(output_S0.shape[1]-1)
        
        loss_corr_fMRInoise = (loss_corr_fMRInoise+abs(torch.mean(torch.abs(corr_output))))/2
        
        """Minimal correlation between BOLD part in brain and non-brain areas"""  
        
        brain_fMRI = y_pred[:,0,:]
        nbrain_fMRI = y_nbrain[:,0,:]
        
        brain_fMRI = brain_fMRI - torch.mean(brain_fMRI,dim=1,keepdim=True)
        brain_fMRI = brain_fMRI / torch.std(brain_fMRI,dim=1,keepdim=True)
        
        nbrain_fMRI = nbrain_fMRI - torch.mean(nbrain_fMRI,dim=1,keepdim=True)
        nbrain_fMRI = nbrain_fMRI / torch.std(nbrain_fMRI,dim=1,keepdim=True)
        
        corr_output = torch.sum(brain_fMRI*nbrain_fMRI,dim=1)/(brain_fMRI.shape[1]-1)
        
        loss_corr_fMRINB = abs(torch.mean(torch.abs(corr_output)))
        
        """BOLD frequencies in range f=[0.008,0.1] Hz"""
        
        fft_BOLD = fft_y_pred[:,0,:]
        
        # Determine which frequencies are included (assuming the rows range from 0Hz to Nyquist)
        f = self.Ny * (np.array(list(range(1, fft_BOLD.shape[1] + 1)))) / (fft_BOLD.shape[1])
        
        # Only include frequencies higher than 0.01Hz and lower than 0.12 Hz
        fincl = np.squeeze(np.array(np.where(np.logical_and(f > 0.008,f < 0.1))))
        
        ALFF = torch.sum(torch.abs(fft_BOLD[:,fincl]),dim=1)
        ftotal = torch.sum(torch.abs(fft_BOLD),dim=1)
        
        fALFF = torch.nan_to_num(ALFF / ftotal,nan=0.0)
        
        loss_fALFF = abs(torch.mean(torch.abs(1 - fALFF)))
        
        """In brain signals: optimal weighted correlation between BOLD part and original echo signals"""
        
        diff_BOLD = 0 
        corr_fMRI = 0 
        loss_nbrain = 0
        for ie in range(len(self.TE)):
            w_bold = (self.TE[ie]/60)*np.exp(-(self.TE[ie]-60)/60)
            w_nbold = (self.TE[ie]/5)*np.exp(-(self.TE[ie]-5)/5)
            num_w_bold = w_bold + w_nbold
            w_bold = w_bold/num_w_bold
            w_nbold = w_nbold/num_w_bold
            
            inputs_gm = y_orig[:,ie,:]          
            BOLD_fMRI = y_pred[:,0,:]
            NBOLD_fMRI = y_pred[:,1,:]
            
            # Difference output - input
            outputs_gm = (BOLD_fMRI * w_bold + NBOLD_fMRI  * w_nbold)
    
            diff_BOLDie = torch.abs(inputs_gm-outputs_gm) 
            diff_BOLD = diff_BOLD + abs(torch.mean(torch.abs(diff_BOLDie))) / len(self.TE)
            
            # Correlation BOLD signal
            BOLD_fMRI = BOLD_fMRI - torch.mean(BOLD_fMRI,dim=1,keepdim=True)
            BOLD_fMRI = BOLD_fMRI / torch.std(BOLD_fMRI,dim=1,keepdim=True)
            
            corr_fMRIie = 1-torch.sum(inputs_gm*BOLD_fMRI,dim=1)/(BOLD_fMRI.shape[1]-1)
            corr_BOLDie = abs(torch.mean(torch.abs(corr_fMRIie))) * w_bold
            
            #Correlation NON-BOLD signal
            NBOLD_fMRI = NBOLD_fMRI - torch.mean(NBOLD_fMRI,dim=1,keepdim=True)
            NBOLD_fMRI = NBOLD_fMRI / torch.std(NBOLD_fMRI,dim=1,keepdim=True)
            
            corr_NfMRIie = 1-torch.sum(inputs_gm*NBOLD_fMRI,dim=1)/(NBOLD_fMRI.shape[1]-1)
            corr_NBOLDie = abs(torch.mean(torch.abs(corr_NfMRIie))) * w_nbold
            
            corr_fMRI = corr_fMRI + (abs(corr_BOLDie) + abs(corr_NBOLDie)) / len(self.TE)
            
        loss_per_echo = (abs(diff_BOLD) + abs(corr_fMRI))/2
        
        return abs(loss_corr_fMRInoise)+abs(loss_corr_fMRINB)+abs(loss_fALFF)+abs(loss_per_echo)