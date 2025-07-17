function AutoSPM1stlevel_fmri

%Script to do the auto 1st level fMRI processing in SPM12
%
%Preparation:
%* Organise the data in BIDS format
%    - datpath
%        -sub-##
%            -ses-00# (if your experiment contains multiple session per subject)
%                -anat: containes the anatomical data (3D T1)
%                -func: containes the fmri data
%            
%* Make sure all data is preprocessed
%    - preproc_anat
%    - preproc_func
    
%* IMPORTANT: !! Look at your preprocessed data before starting any analysis. It makes no sense to lose time in trying to process bad data !!

%% Give path to SPM25

clear 'all'

params.spm_path = '/Users/accurad/Library/Mobile Documents/com~apple~CloudDocs/Matlab/spm25';

%% Give the basic input information of your data

datpath = '/Volumes/LaCie/UZ_Brussel/ASLBOLD_OpenNeuro_FT/IndData'; %'/Volumes/LaCie/UZ_Brussel/ME_fMRI_GE/data';  %'/Volumes/LaCie/UZ_Brussel/ASLBOLD_Manon/data';

sublist = [1:13]; %﻿list with subject id of those to preprocess separated by , (e.g. [1,2,3,4]) or alternatively use sublist = [first_sub:1:last_sub]
params.sub_digits = 2; %if 2 the subject folder is sub-01, if 3 the subject folder is sub-001, ...

nsessions = [1]; %nsessions>0
 
params.task = {'bilateralfingertapping'}; %{'PREcog'}; %{'ME-EmoFaces'}; %text string that is in between task_ and _bold in your fNRI nifiti filename

params.analysisname = '_dune_cbf_glm';
params.modality = 'fasl'; %'fmri' of 'fasl'

params.use_parallel = false; 
params.maxprocesses = 2; %Best not too high to avoid memory problems
params.loadmaxvols = 100; %to reduce memory load, the preprocessing can be split in smaller blocks (default = 100)
params.keeplogs = false;

%% fMRI data parameters
    params.preprocfmridir = 'preproc_dune_wf_asl'; %'preproc_bold_dune'; %'preproc_func_ME-EmoFaces_dune'; %directory with the preprocessed fMRI data
    params.fmri_prefix = 'swcdlavfre'; %'swacdfre'; %'swacdure'; %fMRI file name of form [fmri_prefix 'sub-ii_task-..._' fmri_endfix '.nii']
    
    params.dummytime = 0; %only if the timings in the _events.tsv file should be corrected for dummy scans
        
    %In case of multiple runs in the same session exist
    params.func.mruns = false; %true if run number is in filename
    params.func.runs = [1]; %the index of the runs (in filenames run-(index))
    params.func.use_runs = 'together'; % 'separately' or 'together' (this parameter is ignored if mruns is false)
    %'separately': a separate analysis is done per run
    %'together': all runs are combined in 1 analysis
    
    % For ME-fMRI
    params.func.meepi = true; %true if echo number is in filename
    params.func.echoes = [1:3]; %the index of echoes in ME-fMRI used in the analysis. If meepi=false, echoes=[1]. 

    % For Functional ASL 
    params.whichfile = 'cbf'; %do processing on 'asl' file or on 'cbf' data (default='cbf')
    params.asl.LabelingDuration = 1.525; % in seconds (parameter is ignored if LabelingDuration is in json file)
    params.asl.PostLabelDelay = 1.525; % in seconds (parameter is ignored if PostLabelDelay is in json file)

%% SPM first level analysis parameters
    params.analysis_type = 'GLM'; % 'GLM' or 'within_subject' (default='GLM'), 'within_subject' only possible for 'fasl'

    params.confounds_prefix = 'rp_e'; %confounds file of form [confounds_prefix 'sub-ii_task-... .txt']
    params.add_parametricModulation = false; %use the weights in events.tsv for parametric modultion
    params.add_regressors = false; %if data not denoised set true otherwhise false 
    params.add_derivatives = true; %add temmperal and dispertion derivatives to the GLM (default=false)
    params.optimize_HRF = false; %Optimize HRF parameters (peak time and duration) to the data using the TEDM toolbox (only possible for BOLD)
    params.use_ownmask = true;
    params.model_serial_correlations = 'AR(1)'; %'AR(1) for fmri, 'none' for fasl
    params.hpf = 128; %default 128 but changed to tr*(nvol-1) if already filtered (f in prefix)

%% SPM results analysis
    %Save SPM results per ccontrast as thresholded map, binary mask, n-aray map (n=cluster number), 
    %csv file, pdf file
    params.save_spm_results = true;
    params.threshold_correction = 'none'; %'none' or 'FWE' (default='none')
    params.pthreshold = 0.001; % signifiance p-threshold (default=0.001) 
    params.kthreshold = 0; %ccluster extend threshold (default=0)
    params.save_thresholded_map = true;
    params.save_binary_mask = true;
    params.save_naray = false;
    params.save_csv_file = false;
    params.save_pdf_file = false;
    params.save_tiff_file = true;

%% Define the contrasts
    %contrast(i) is structure with fields
    %   conditions: conditions to compare
    %   vector: contrast weight vector
    %e.g A contrast between conditions is given as
    %   contrast(i).conditions={'condition 1','condition 2'};
    %   contrast(i).vector=[1 -1];

    %params.contrast(1).conditions = {'Go_Left','Go_Right','NoGo_Left','NoGo_Right'};
    %params.contrast(1).vector = [1,1,-1,-1];

    %params.contrast(2).conditions = {'Go_Left','Go_Right','NoGo_Left','NoGo_Right'};
    %params.contrast(2).vector = [-1,-1,1,1];

    params.contrast(1).conditions = {'Finger'};
    params.contrast(1).vector = [1];

    params.contrast(2).conditions = {'Finger'};
    params.contrast(2).vector = [-1];

    %params.contrast(1).conditions = {'episodic','semantic'};
    %params.contrast(1).vector = [1,-1];
    
    %params.contrast(2).conditions = {'episodic','semantic'};
    %params.contrast(2).vector = [-1,1];
    
    %params.contrast(1).conditions = {'sad','neutral'};
    %params.contrast(1).vector = [1,-1];
    
    %params.contrast(2).conditions = {'sad','neutral'};
    %params.contrast(2).vector = [-1,1];
    
    %params.contrast(3).conditions = {'happy','neutral'};
    %params.contrast(3).vector = [1,-1];
    
    %params.contrast(4).conditions = {'happy','neutral'};
    %params.contrast(4).vector = [-1,1];
    
    %params.contrast(5).conditions = {'sad','happy'};
    %params.contrast(5).vector = [1,-1];
    
    %params.contrast(6).conditions = {'sad','happy'};
    %params.contrast(6).vector = [-1,1];

%% BE CAREFUL WITH CHANGING THE CODE BELOW THIS LINE !!
%--------------------------------------------------------------------------------

restoredefaultpath

[params.my_spmbatch_path,~,~] = fileparts(mfilename('fullpath'));

if exist(params.spm_path,'dir'), addpath(genpath(params.spm_path)); end
if exist(params.my_spmbatch_path,'dir'), addpath(genpath(params.my_spmbatch_path)); end

old_spm_read_vols_file=fullfile(spm('Dir'),'spm_read_vols.m');
new_spm_read_vols_file=fullfile(spm('Dir'),'old_spm_read_vols.m');

if isfile(old_spm_read_vols_file), movefile(old_spm_read_vols_file,new_spm_read_vols_file); end
  
fprintf('Start with processing the data\n')

warnstate = warning;
warning off;

% User interface.
SPMid                 = spm('FnBanner',mfilename,'2.10');
[Finter,Graf,CmdLine] = spm('FnUIsetup','Preproces SPM');

spm('defaults', 'FMRI');

curdir = pwd;

my_spmbatch_start_fmriprocessing(sublist,nsessions,datpath,params);

spm_figure('close',allchild(0));

cd(curdir)

if isfile(new_spm_read_vols_file), movefile(new_spm_read_vols_file,old_spm_read_vols_file); end

fprintf('\nDone with processing the data\n')

clear 'all'

end