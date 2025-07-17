function AutoDCM2nii

%Script to do the auto convertion from DCM file into nii files
%
%* The nii data will be saved in BIDS format
%    - datpath
%        -sub-##
%            -ses-00# (if your experiment contains multiple session per subject)
%                -anat: containes the anatomical data (3D T1)
%                   Files: sub-##_T1w.nii and sub-##_T1w.json
%                -func: containes the fmri data
%                   Files: sub-##_task-..._bold.nii and sub-##_task-..._bold.json
%                -fmap: containnes the gradient pololarity (blip-up/down) filpt data or the fieldmap scans
%                   Files in case of inverted gradient polarity: sub-##_dir-pi_epi.nii and sub-##_dir-pi_epi.json
%                   Files in case of fieldmap scans: (image 1 in file is amplitude, image 2 in file is phase)
%                          sub-##_fmap_echo-1.nii and sub-##_fmap_echo-1.json
%                          sub-##_fmap_echo-2.nii and sub-##_fmap_echo-2.json
    
%* IMPORTANT: !! Look at your data before starting any (pre)processing. Losing time in trying to process bad data makes no sense !!

%Script written by dr. Peter Van Schuerbeek (Radiology UZ Brussel)

%% Give path to SPM25

params.spm_path = '/Users/accurad/Library/Mobile Documents/com~apple~CloudDocs/Matlab/spm25';

%% Give the basic input information of your data

params.datpath = '/Volumes/LaCie/UZ_Brussel/ASLBOLD_Manon/RawData';

sublist = [1];%list with subject id of those to preprocess separated by , (e.g. [1,2,3,4]) or alternatively use sublist = [first_sub:1:last_sub]
params.sub_digits = 2; %if 2 the result will be sub-01, if 3 the result will be sub-001

%Add per sequence to convert an extra ssequence object to the mri_data structure as (folder,seqtype,name,task,[session],add_run,add_echo,add_acq)
% folder: substructure starting from sub-## containing the dicom files
% acqtype: (anat, func, fmap, dti, perf) will be used as foldername to write the converted nifti files
% seqtype: type of sequence to make the name (T1w, pepolar, fmri (bold), aslbold, pcasl)
% task: for fMRI task-string (what comes in the name as _task-...) / if seqtype not fmri, set ''
% session: scan session
% run: run number
% add_acq: (True or False) add sequence name to nifti file name (_acq-...)  
% add_dir: (True or False) add phase encoding direction to nifti file name (_dir-...)
% add_run: (True or False) add run numer tag to nifti file name (_run-#)
% add_echo: (True or False) add echo numer tag to nifti file name (_echo-#)

%HOW TO:
% 1. To add a sequence, copy a params.mridata(ii) blocck and set the parameters
% 2. params.mridata(ii) ii should start at 1
% 3. To remove a params.mridata(ii) block, by deleting the block or by
% deactivating the lines by addinng a % at the beginning

%% Example anatomical T1w scan
%params.mridata(1).folder = 'ses-001/VBM';
%params.mridata(1).acqtype = 'anat';
%params.mridata(1).seqtype = 'T1w';
%params.mridata(1).task = '';
%params.mridata(1).session = 1;
%params.mridata(1).run = 1;
%params.mridata(1).add_acq = false;
%params.mridata(1).add_dir = false;
%params.mridata(1).add_run = false;
%params.mridata(1).add_echo = false;

%% Example pepolar fmap scan
params.mridata(1).folder = 'ses-001/PP_ASLBOLD';
params.mridata(1).acqtype = 'fmap';
params.mridata(1).seqtype = 'pepolar';
params.mridata(1).task = '';
params.mridata(1).session = 1;
params.mridata(1).run = 1;
params.mridata(1).add_acq = false;
params.mridata(1).add_dir = true;
params.mridata(1).add_run = false;
params.mridata(1).add_echo = true;

%% Example fMRI scan
%params.mridata(1).folder = 'dcm/SE-fMRI_EFT';
%params.mridata(1).acqtype = 'func';
%params.mridata(1).seqtype = 'fmri';
%params.mridata(1).task = 'SE-EFT';
%params.mridata(1).session = 1;
%params.mridata(1).run = 1;
%params.mridata(1).add_acq = false;
%params.mridata(1).add_dir = true;
%params.mridata(1).add_run = false;
%params.mridata(1).add_echo = true;

%params.mridata(2).folder = 'dcm/SE-fMRI_EmoFaces';
%params.mridata(2).acqtype = 'func';
%params.mridata(2).seqtype = 'fmri';
%params.mridata(2).task = 'SE-EmoFaces';
%params.mridata(2).session = 1;
%params.mridata(2).run = 1;
%params.mridata(2).add_acq = false;
%params.mridata(2).add_dir = true;
%params.mridata(2).add_run = false;
%params.mridata(2).add_echo = true;

%% Example ASLBOLD scan
%params.mridata(1).folder = 'ses-001/STROOP1';
%params.mridata(1).acqtype = 'func';
%params.mridata(1).seqtype = 'aslbold';
%params.mridata(1).task = 'stroop';
%params.mridata(1).session = 1;
%params.mridata(1).run = 1;
%params.mridata(1).add_acq = false;
%params.mridata(1).add_dir = true;
%params.mridata(1).add_run = true;
%params.mridata(1).add_echo = true;

%params.mridata(4).folder = 'ses-001/PREcog';
%params.mridata(4).acqtype = 'func';
%params.mridata(4).seqtype = 'aslbold';
%params.mridata(4).task = 'PREcog';
%params.mridata(4).session = 1;
%params.mridata(4).run = 1;
%params.mridata(4).add_acq = false;
%params.mridata(4).add_dir = true;
%params.mridata(4).add_run = true;
%params.mridata(4).add_echo = true;

params.use_parallel = true; %only possible when parallel toolbox is installed

%Be careful with changing the code below this line!
%--------------------------------------------------------------------------------
restoredefaultpath

[params.my_spmbatch_path,~,~] = fileparts(mfilename('fullpath'));

if exist(params.spm_path,'dir'), addpath(genpath(params.spm_path)); end
if exist(params.my_spmbatch_path,'dir'), addpath(genpath(params.my_spmbatch_path)); end

old_spm_read_vols_file=fullfile(spm('Dir'),'spm_read_vols.m');
new_spm_read_vols_file=fullfile(spm('Dir'),'old_spm_read_vols.m');

if isfile(old_spm_read_vols_file), movefile(old_spm_read_vols_file,new_spm_read_vols_file); end

fprintf('Start with converting the data\n')

curdir = pwd;

warnstate = warning;
warning off;

% User interface.
SPMid                 = spm('FnBanner',mfilename,'2.10');
[Finter,Graf,CmdLine] = spm('FnUIsetup','Preproces SPM');

spm('defaults', 'FMRI');

curdir = pwd;

my_spmbatch_noparalleldcm2nii(sublist,params)

spm_figure('close',allchild(0));

cd(curdir)

if isfile(new_spm_read_vols_file), movefile(new_spm_read_vols_file,old_spm_read_vols_file); end

clear 'all'

fprintf('\nDone with converting the data\n')

end