function AutoSPMpreprocessing_aslbold

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
%                -func: containes the fmri data (BOLD only or ASL-BOLD)
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

datpath = '/Volumes/LaCie/UZ_Brussel/ASLBOLD_OpenNeuro_FT/IndData';

sublist = [1:13];%list with subject id of those to preprocess separated by , (e.g. [1,2,3,4]) or alternatively use sublist = [first_sub:1:last_sub]
params.sub_digits = 2; %if 2 the subject folder is sub-01, if 3 the subject folder is sub-001, ...

nsessions = [1]; %nsessions>0

params.func_save_folder = 'preproc_dune_wf_bold'; %name of the folder to save the preprocessed bold data
params.perf_save_folder = 'preproc_dune_wf_asl'; %name of the folder to save the preprocessed asl data

task ={'bilateralfingertapping'};

%% In case of multiple runs in the same session exist
params.func.mruns = false; %true if run number is in filename
params.func.runs = [1]; %the index of the runs (in filenames run-(index))

%% Parallel processing and memory reduction
params.use_parallel = true; %(default=false)
params.maxprocesses = 2; %Best not too high to avoid memory problems %(default=2)
params.loadmaxvols = 100; %to reduce memory load, the preprocessing can be split in smaller blocks (default = 100)
params.keeplogs = false; %(default=false)

%% Save intermediate results needed?
params.save_intermediate_results = true; %clean up the directory by deleting unnecessary files generated during the processing (default = false)

%% Which analyses to do
params.preprocess_anatomical = false;  %(default=true)  
params.preprocess_functional = true; %(default=true)
params.preprocess_asl = true; %(default=true)

%% FMRI parameters
params.func.meepi = true; %true if echo number is in filename (default=true)
params.func.echoes = [1:4]; %the index of echoes in ME-fMRI used in the analysis. If meepi=false, echoes=[1]. 

params.func.dummytime = 0; %time in seconds (default=2*TR)

params.func.pepolar = false; %true if fmap scan exist otherwise false (default=true)

%% ASL Parameters
params.asl.splitaslbold = 'dune'; %'meica' or 'dune' (default='meica') 
%'meica': after filtering, ME-ICA (tedana based)
%'dune': experimental splitting method

%---------------------------------------------------------------------------------------
%%% BE CAREFUL WITH CHANGING THE CODE BELOW THIS LINE !!
%---------------------------------------------------------------------------------------

%% Preprocessing anatomical data

    % Normalization
    params.anat.do_normalization = true;
    params.anat.normvox = [2.0 2.0 2.0]; %(default=[2.0 2.0 2.0]) Same as for fMRI!!

    % Segmentation using CAT12
    params.anat.do_segmentation = false; %(default=true)
    params.anat.roi_atlas = false; %(default=false)
    
%% Preprocessing ASL data

    %ASL data
    params.asl.isM0scan = 'last'; %The M0 image is by defaullt the last volume @GE (set 'last')
    params.asl.LabelingDuration = 1.525; % in seconds (parameter is ignored if LabelingDuration is in json file)
    params.asl.PostLabelDelay = 1.525; % in seconds (parameter is ignored if PostLabelDelay is in json file)

    params.asl.GMWM = 'anat'; %wich data used for segmentation maps 'anat', 'M0asl' (default='anat')
    
%% Preprocessing functional (the order of the parameters represents the fixed order of the steps done)

    % Realignnment (motion correction)
    params.func.do_realignment = true; %(default=true)
    
    % Detect and correct bad volumes
    params.func.do_ArtRepair = true; %(default=true)

    %Denoising before echo combination and normalization ussing the
    %parameters from params.denoise
    params.func.denoise = true; %(default=false)

    params.func.combination = 'T2star_weighted'; %only used for ME-EPI (default=T2star_weighted)
    %none: all echoes are preprocessed separatly
    %average: The combination is the average of the multiple echo images
    %TE_weighted: The combination is done wi=TEi or 
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

%% Start the analysis
    
restoredefaultpath

[params.my_spmbatch_path,~,~] = fileparts(mfilename('fullpath'));

if exist(params.spm_path,'dir'), addpath(genpath(params.spm_path)); end
if exist(params.GroupICAT_path,'dir'), addpath(genpath(params.GroupICAT_path)); end
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

my_spmbatch_start_aslboldpreprocessing(sublist,nsessions,task,datpath,params)

spm_figure('close',allchild(0));

cd(curdir)

if isfile(new_spm_read_vols_file), movefile(new_spm_read_vols_file,old_spm_read_vols_file); end

fprintf('\nDone\n')

clear 'all'

end