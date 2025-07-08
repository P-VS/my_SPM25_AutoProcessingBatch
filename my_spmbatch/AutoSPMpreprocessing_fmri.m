function AutoSPMpreprocessing_fmri

%Script to do the auto preprocessing in SPM12
%
%Preparation:
% Convert the DICOM files into nifti using dcm2niix in MROCroGL
% For the anatomical scans, set 'Crop 3D Images' on
%
%* Organise the data in BIDS format
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

%% Give path to SPM25 and GroupICA

params.spm_path = '/Users/accurad/Library/Mobile Documents/com~apple~CloudDocs/Matlab/spm25';
params.GroupICAT_path = '/Users/accurad/Library/Mobile Documents/com~apple~CloudDocs/Matlab/GroupICATv40c';

%% Give the basic input information of your data

datpath = '/Volumes/LaCie/UZ_Brussel/ME_fMRI_GE/data';

sublist = [1];%list with subject id of those to preprocess separated by , (e.g. [1,2,3,4]) or alternatively use sublist = [first_sub:1:last_sub]
params.sub_digits = 2; %if 2 the subject folder is sub-01, if 3 the subject folder is sub-001, ...

nsessions = [1]; %nsessions>0

params.func_save_folder = 'preproc_func_SPM25'; %name of the folder to save the preprocessed bold data

task ={'ME-EFT'};

%In case of multiple runs in the same session exist
params.func.mruns = false; %true if run number is in filename
params.func.runs = [1]; %the index of the runs (in filenames run-(index))

% For ME-fMRI
params.func.meepi = true; %true if echo number is in filename
params.func.echoes = [1]; %the index of echoes in ME-fMRI used in the analysis. If meepi=false, echoes=[1]. 

params.use_parallel = false; 
params.maxprocesses = 2; %Best not too high to avoid memory problems
params.loadmaxvols = 100; %to reduce memory load, the preprocessing can be split in smaller blocks (default = 1000)
params.keeplogs = false;

params.save_intermediate_results = true; %clean up the directory by deleting unnecessary files generated during the processing (default = false)

params.preprocess_anatomical = false;
params.preprocess_functional = true;

% Geometric correction if fmap is available
params.func.pepolar = true; %(default=true)
   
%---------------------------------------------------------------------------------------
% BE CAREFUL WITH CHANGING THE CODE BELOW THIS LINE !!
%---------------------------------------------------------------------------------------

%% Preprocessing anatomical data

    % Normalization
    params.anat.do_normalization = true; %(default=true)
    params.anat.normvox = [2.0 2.0 2.0]; %(default=[2.0 2.0 2.0]) Same as for fMRI!!

    % Segmentation ussing CAT12
    params.anat.do_segmentation = false; %(default=true)
    params.anat.roi_atlas = false; %(default=false)
    
%% Preprocessing functional data (the order of the parameters represents the fixed order of the steps done)

    % Remove the dummy scans n_dummy_scans = floor(dummytime/TR)
    params.func.dummytime = 8; %time in seconds
    
    % Realignnment (motion correction)
    params.func.do_realignment = true; %(default=true)

    % Detect and correct bad volumes
    params.func.do_ArtRepair = true; %(default=true)

    %Denoising before echo combination and normalization ussing the
    %parameters from params.denoise
    params.func.denoise = true; %(default=true)
 
    params.func.combination = 'T2star_weighted'; %only used for ME-EPI (default=T2star_weighted)
    %none: all echoes are preprocessed separatly
    %average: The combination is the average of the multiple echo images
    %TE_weighted: The combination is done wi=TEi
    %T2star_weighted: dynamic T2* weighted combination (loglinear fit from tedana)
    %dyn_T2star: dynamic T2* mapping (not adviced due to the propagation of noise leading to spikes in the T2* maps)
    %see Heunis et al. 2021. The effects of multi-echo fMRI combination and rapid T2*-mapping on offline and real-time BOLD sensitivity. NeuroImage 238, 118244
         
    % Slice time correction
    params.func.do_slicetime = true; %(default=true)
      
    % Normalization
    params.func.do_normalization = true; %(default=true)
    params.func.normvox = [2.0 2.0 2.0]; %(default=[2.0 2.0 2.0])
     
    % Smoothing
    params.func.do_smoothing = true; %(default=true)
    params.func.smoothfwhm = 6; %(default=6)

%% Denoising (after normalization)
    %if do_denoising = true and preprocess_functional = false -> only denoising of the normalised data
    %Denoising will not run if it was already done before

    params.do_denoising = false; %(default=false)

    % Extend motion regressors with derivatives and squared regressors
    params.denoise.do_mot_derivatives = true; %derivatives+squares (24 regressors) (default=true)

    % Band-pass filtering
    params.denoise.do_bpfilter = true; %(default=true)
    params.denoise.bpfilter = [0.008 Inf]; %no highpass filter is first 0, no lowpass filter is last Inf, default=[0.008 Inf]
    params.denoise.polort = 1; %order of the polynomial function used to remove the signal trend (0: only mean, 1: linear trend, 2: quadratic trend, default=2)

    % aCompCor
    params.denoise.do_aCompCor = true; %(default=true)
    params.denoise.Ncomponents = 5; %if in range [0 1] then the number of aCompCor components is equal to the number of components that explain the specified percentage of variation in the signal (default=5)

    % ICA-AROMA
    params.denoise.do_ICA_AROMA = true; %(default=true)

    % Noise regression / remove ICA-AROMA noise components
    params.denoise.do_noiseregression = false; %(default=true)

    % Prepare data for DUNE denoising in python (WIP)
    params.denoise.do_DUNE = false; %(default=false)
    params.denoise.DUNE_part = 'bold'; %'bold' or 'nonbold' (default='bold')
   
%% Start preprocessing

restoredefaultpath

[params.my_spmbatch_path,~,~] = fileparts(mfilename('fullpath'));

if exist(params.GroupICAT_path,'dir'), addpath(genpath(params.GroupICAT_path)); end
if exist(params.spm_path,'dir'), addpath(genpath(params.spm_path)); end
if exist(params.my_spmbatch_path,'dir'), addpath(genpath(params.my_spmbatch_path)); end

old_spm_read_vols_file=fullfile(spm('Dir'),'spm_read_vols.m');
new_spm_read_vols_file=fullfile(spm('Dir'),'old_spm_read_vols.m');

if isfile(old_spm_read_vols_file), movefile(old_spm_read_vols_file,new_spm_read_vols_file); end
  
fprintf('Start with preprocessing \n')

curdir = pwd;

warnstate = warning;
warning off;

% User interface.
SPMid                 = spm('FnBanner',mfilename,'2.10');
[Finter,Graf,CmdLine] = spm('FnUIsetup','Preproces SPM');

spm('defaults', 'FMRI');

my_spmbatch_start_fmripreprocessing(sublist,nsessions,task,datpath,params)

spm_figure('close',allchild(0));

cd(curdir)

if isfile(new_spm_read_vols_file), movefile(new_spm_read_vols_file,old_spm_read_vols_file); end

if exist(fullfile(datpath,'derivatives')), rmdir(fullfile(datpath,'derivatives'),'s'); end

fprintf('\nDone\n')

clear all

end