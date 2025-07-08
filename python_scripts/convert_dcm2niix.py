#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 13 13:14:04 2021

@author: dr. Peter Van Schuerbeek
"""

"""
Converting DICOM into .nii data

Organise the data in BIDS format
    - datpath
        -sub-##
            -ses-00# (if your experiment contains multiple session per subject)
                -anat: containes the anatomical data (3D T1)
                   Files: sub-##_T1w.nii and sub-##_T1w.json
                -func: containes the fmri data
                   Files: sub-##_task-..._bold.nii and sub-##_task-..._bold.json
                -fmap: containnes the gradient pololarity (blip-up/down) filpt data or the fieldmap scans
                   Files in case of inverted gradient polarity: sub-##_dir-##_epi.nii and sub-##_dir-pi_epi.json
                   Files in case of fieldmap scans: sub-##_phase1.nii, sub-##_phase2.nii, sub-##_magnitude1.nii and sub-##_magnitude2.nii
                -perf: ASL data based on GE PCASL sequence which gives only the m0, deltam and CBF images
                    Files: sub-##_asl_m0scan.nii
                    Files: sub-##_asl_deltam.nii
                    Files: sub-##_asl_cbf.nii
"""

import warnings
import sys
import os

if not sys.warnoptions:
    warnings.simplefilter("ignore")

import shutil
import json

from nipype import Workflow, Node, IdentityInterface

from nipype.interfaces.io import SelectFiles
from nipype.interfaces.utility import Function
from nipype.interfaces.dcm2nii import Dcm2niix

def set_preprocessing_parameters():
    
    """
    Give the basic input information of your data
    """
    pp_params = {}
    
    pp_params['datpath'] = '/Volumes/LaCie/UZ_Brussel/ASLBOLD_Manon/RawData'  #No spaties in the path name
    
    first_sub = 1
    last_sub = 1
    pp_params['sublist'] = [2] #list with subject id of those to preprocess separated by , (e.g. [1,2,3,4]) or alternatively use sublist = list(range(first_sub,last_sub+1))
    pp_params['sub_digits'] = 2 #if 2 the result will be sub-01, if 3 the result will be sub-001
    
    """
    Add per sequence to convert an extra ssequence object to the mri_data structure as (folder,seqtype,name,task,[session],add_run,add_echo,add_acq)
    0.folder: substructure starting from sub-## containing the dicom files
    1.acqtype: (anat, func, fmap, dti, perf) will be used as foldername to write the converted nifti files
    2.seqtype: type of sequence to make the name (T1w, fieldmap, pepolar, fmri, asl)
    3.task: for fMRI task-string (what comes in the name as _task-...) / if seqtype not fmri, set ''
    4.session: scan session
    5.run: run number
    6.add_acq: (True or False) add sequence name to nifti file name (_acq-...)  
    7.add_dir: (True or False) add phase encoding direction to nifti file name (_dir-...)
    8.add_run: (True or False) add run numer tag to nifti file name (_run-#)
    9.add_echo: (True or False) add echo numer tag to nifti file name (_echo-#)
    
    example
    T1w:      scan_1 = ('dicom/anat_t1','anat','T1w','',[1],[1],False,False,False,False)
    fieldamp: scan_1 = ('dicom/fieldmap','fmap','fieldmap','',[1],[1],False,False,False,False) !!based on 2 GE scans TE_1 and TE_2!!
    pepolar:  scan_1 = ('dicom/fmri_pp','fmap','pepolar','',[1],[1],False,True,False,False)
    fMRI:     scan_1 = ('dicom/fmri_task','func','fmri','rest',[1],[1],False,True,True,True)
    ASLBOLD:  scan_1 = ('dicom/aslbold','func','aslbold','rest',[1],[1],False,True,True,True)
    pcasl:    scan_1 = ('dicom/pcasl','perf','asl','',[1],[1],False,False,False,False) !!GE PCASL based!!
    """
    
    pp_params['mri_data'] = []
    
    scan_1 = ('ses-001/VBM','anat','T1w','',[1],[1],False,False,False,False)
    scan_2 = ('ses-001/STROOP1','func','aslbold','test',[1],[1],False,True,True,True)
    
    pp_params['mri_data'].append([scan_1,scan_2])
    
    return pp_params


"""
BE CAREFUL WITH CHANGING THE CODE BELOW THIS LINE !!
---------------------------------------------------------------------------------------
"""

def load_d2n_params(subdir,folder,acqtype,seqtype,task,session,run):
    
    import os
    
    split_path = subdir.split(os.sep)
    substring = split_path[len(split_path)-1]
    
    sesstring = 'ses-00'+str(session)
    
    source_dir = os.path.join(subdir,folder)
    output_dir = os.path.join(subdir,sesstring,acqtype)
    
    out_filename = substring+'_'+sesstring
    
    if 'fmri' in seqtype:
        out_filename = out_filename+'_task-'+task
        
    out_filename = out_filename+'_acq-'+'%d'
    out_filename = out_filename+'_run-'+str(run)
    out_filename = out_filename+'_echo-'+'%e'
    
    crop = False
    ignore_deriv = False
    
    if 'anat' in acqtype:
        crop = True
        ignore_deriv = True
    
    return source_dir, crop, ignore_deriv, out_filename, output_dir

"""
---------------------------------------------------------------------------------------
"""
def check_zipfile(source_dir):
    
    import os
    import shutil
    
    ldir = os.listdir(source_dir)
    
    for ifile in ldir:
        if os.path.isfile(os.path.join(source_dir,ifile)):
            if '.zip' in ifile and not '._' in ifile:
                shutil.unpack_archive(os.path.join(source_dir,ifile),source_dir)
                       
    return source_dir
    
    
"""
---------------------------------------------------------------------------------------
"""
def rename_file(in_files,seqtype,add_acq,add_dir,add_run,add_echo):
    
    import os
    import json
    
    if 'pepolar' in seqtype: add_dir=True
    
    if len(in_files[1])==1: in_files = [in_files]
   
    for ifile in in_files:
        ifile = ifile.split('.nii')[0]
        acq_split = ifile.split('_acq-')
        run_split = acq_split[1].split('_run-')
        echo_split = run_split[1].split('_echo-')  
        ph_split = echo_split[1].split('_ph')
        new_file = acq_split[0]
        
        if add_acq: new_file = new_file+'_acq-'+run_split[0]
        
        if add_dir:
            with open(ifile+'.json','r') as f: jsondat=json.load(f)
            f.close()
            
            if "PhaseEncodingDirection" in jsondat:
                pedir = jsondat["PhaseEncodingDirection"]
            else: 
                pedir = ''
                pestring = 'NF' #Not found
            
            if 'i' in pedir: pestring = 'LR'
            if 'j' in pedir: pestring = 'PA'
            if 'k' in pedir: pestring = 'FH'
            
            if '-' in pedir: pestring = pestring[::-1]
            
            new_file = new_file+'_dir-'+pestring
  
        if add_run: new_file = new_file+'_run-'+echo_split[0]
      
        if not add_echo:
            test_string = acq_split[0]+'_acq-'+run_split[0]+'_run-'+echo_split[0]+'_echo-2'
            
            test = [i for i in in_files if test_string in i]
            if len(test)>0: add_echo=True
            
        if add_echo: new_file = new_file+'_echo-'+ph_split[0]
        
        if 'T1w' in seqtype: 
            new_file = new_file+'_T1w'
            
            if os.path.isfile(ifile+'_Crop_1.nii'):
                os.rename(ifile+'_Crop_1.nii',new_file+'_Crop_1.nii')
        
        if 'fieldmap' in seqtype: #In some older studies done on our GE scanner: 2 GE scans at different TEs (no longer used)
            if '_ph' in echo_split[1]: #Based on sequence names
                if 'TE_1' in acq_split[1]: new_file = new_file+'_phase1'
                if 'TE_2' in acq_split[1]: new_file = new_file+'_phase2'
            else:
                if 'TE_1' in acq_split[1]: new_file = new_file+'_magnitude1'
                if 'TE_2' in acq_split[1]: new_file = new_file+'_magnitude2'
        
        if 'asl' in seqtype:  
            with open(ifile+'.json','r') as f: jsondat=json.load(f)
            f.close()
                        
            ImType = jsondat["ImageType"]
            
            if 'ORIGINAL' in ImType[0]:
                ftype = '_m0scan'
            elif 'DERIVED' in ImType[0]:
                if 'PRIMARY' in ImType[1]:
                    ftype = '_deltam'
                elif 'SECONDARY' in ImType[1]:
                    ftype = '_cbf'
                    
            new_file = new_file+ftype
            
        if 'pepolar' in seqtype: new_file = new_file+'_epi'
        if 'fmri' in seqtype: new_file = new_file+'_bold' 
        if 'aslbold' in seqtype: new_file = new_file+'_aslbold' 
            
        os.rename(ifile+'.nii',new_file+'.nii')
        os.rename(ifile+'.json',new_file+'.json')
        
    
    return '' #out_files

"""
---------------------------------------------------------------------------------------
"""

def main():
    
    pp_params = set_preprocessing_parameters()
    
    datpath = pp_params['datpath']

    sublist = pp_params['sublist']
    sub_digits = pp_params['sub_digits']
    
    mri_data = pp_params['mri_data']
    
    substringslist = list()
    
    for i in sublist:

        if i<10:
            if sub_digits==2: substring = 'sub-0'+str(i)
            if sub_digits==3: substring = 'sub-00'+str(i)
        elif i<100:
            if sub_digits==2: substring = 'sub-'+str(i)
            if sub_digits==3: substring = 'sub-0'+str(i)
        else:
            substring = 'sub-'+str(i)
            
        substringslist.append(substring)
        
        for k in range(0,len(mri_data[0])):
            output_dir = os.path.join(datpath,substring,'ses-00'+str(mri_data[0][k][4][0]))
            if not os.path.isdir(output_dir): os.mkdir(output_dir) 
            
            output_dir = os.path.join(datpath,substring,'ses-00'+str(mri_data[0][k][4][0]),mri_data[0][k][1])
            if not os.path.isdir(output_dir): os.mkdir(output_dir) 
            
    
    templates = {}
    
    templates['subdir'] = os.path.join(datpath,'{substring}')
    
    """
    Create a preprocessing workflow and select input files
    """
    print('Make workflow step: initiate')
    
    infosource = Node(IdentityInterface(fields=['substring']),name='infosource')

    infosource.iterables = [('substring', substringslist)]

    dcm2niiwf = Workflow(base_dir=datpath,name='dcm2niiwf')
    
    selectfiles = Node(SelectFiles(templates,base_directory=datpath),name="selectfiles")
    
    dcm2niiwf.connect(infosource, 'substring', selectfiles, 'substring')
    
    for k in range(0,len(mri_data[0])):
        load_par_node = Node(interface=Function(input_names=['subdir','folder','acqtype','seqtype','task','session','run'],
                                                output_names=['source_dir','crop','ignore_deriv','out_filename','output_dir'],
                                                function=load_d2n_params),name='load_par'+str(k))
        
        load_par_node.inputs.folder = mri_data[0][k][0]
        load_par_node.inputs.acqtype = mri_data[0][k][1]
        load_par_node.inputs.seqtype = mri_data[0][k][2]
        load_par_node.inputs.task = mri_data[0][k][3]
        load_par_node.inputs.session = mri_data[0][k][4][0]
        load_par_node.inputs.run = mri_data[0][k][5][0]
               
        dcm2niiwf.connect([(selectfiles,load_par_node,[('subdir','subdir')])])
        
        checkzip_node = Node(interface=Function(input_names=['source_dir'],
                                                output_names=['source_dir'],
                                                function=check_zipfile),name='check_zip'+str(k))
        
        dcm2niiwf.connect([(load_par_node,checkzip_node,[('source_dir','source_dir')])])
        
        d2nii_node = Node(Dcm2niix(compress='n', anon_bids=True, bids_format=True), name='d2nii'+str(k))
        
        dcm2niiwf.connect([(load_par_node,d2nii_node,[('crop','crop'),
                                                      ('ignore_deriv','ignore_deriv'),
                                                      ('out_filename','out_filename'),
                                                      ('output_dir','output_dir')
                                                      ]),
                           (checkzip_node,d2nii_node,[('source_dir','source_dir')])
                           ])
        
        renamefile_node = Node(interface=Function(input_names=['in_files','seqtype','add_acq','add_dir','add_run','add_echo'],
                                                  output_names=['out_files'],
                                                  function=rename_file),name='rename_file'+str(k))
        
        renamefile_node.inputs.seqtype = mri_data[0][k][2]
        renamefile_node.inputs.add_acq = mri_data[0][k][6]
        renamefile_node.inputs.add_dir = mri_data[0][k][7]
        renamefile_node.inputs.add_run = mri_data[0][k][8]
        renamefile_node.inputs.add_echo = mri_data[0][k][9]
        
        dcm2niiwf.connect([(d2nii_node,renamefile_node,[('converted_files','in_files')])])
    
    """
    Run the workflow
    """

    print('')
    print('Start dcm2nii sub '+str(i))
    print('')
    dcm2niiwf.run()
    #dcm2niiwf.run(plugin='MultiProc')
    print('Done dcm2nii sub '+str(i))
    print('')

    shutil.rmtree(os.path.join(datpath,'dcm2niiwf'), ignore_errors=True)
 
if __name__ == '__main__':
    main()