function AutoSPMpreprocessing_asl

%Script to do the auto ASL preprocessing in SPM12
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
%                -perf: containes the ASL data
%                   Files: sub-##_asl_m0scan.nii and sub-##_asl_m0scan.json
%                          sub-##_asl_deltam.nii and sub-##_asl_deltam.json
%                          sub-##_asl_cbf.nii and sub-##_asl_cbf.json
%    
%* IMPORTANT: !! Look at your data before starting any (pre)processing. Losing time in trying to process bad data makes no sense !!

%Script written by dr. Peter Van Schuerbeek (Radiology UZ Brussel)

%% Give path to SPM25

params.spm_path = '/Users/accurad/Library/Mobile Documents/com~apple~CloudDocs/Matlab/spm25';

%% Give the basic input information of your data

datpath = '/Volumes/LaCie/UZ_Brussel/rTMS-fMRI_Interleaved/Data';

sublist = [3];%list with subject id of those to preprocess separated by , (e.g. [1,2,3,4]) or alternatively use sublist = [first_sub:1:last_sub]
nsessions = [2]; %nsessions>0

params.save_folder = 'preproc_aslge_test';

params.use_parallel = false; 
params.maxprocesses = 4; %Best not too high to avoid memory problems
params.keeplogs = false;

params.save_intermediate_results = false; 

params.reorient = true; % align data with MNI template to improve normalization

%% Preprocessing anatomical data

    params.preprocess_anatomical = false;

    % Normalization
    params.anat.do_normalization = true;
    params.anat.normvox = [2 2 2];

    % Segmentation ussing CAT12
    params.anat.do_segmentation = false;
    params.anat.roi_atlas = false;
    
%% Preprocessing the PCASL data (based on GE PCASL scans)

    params.preprocess_pcasl = true;

    %In case of multiple runs in the same session exist
    params.asl.mruns = false; %true if run number is in filename
    params.asl.runs = [1]; %the index of the runs (in filenames run-(index))
    
    % In case data already processsed at the GE scanner (PCASL scan)
    params.asl.do_cbfmapping = true;
    params.asl.T1correctionM0 = 'tisssue_maps'; %T1 correction of M0scan: 'tisssue_maps', 'T1_map', 'average_GM', 'average_WM'

    % Normalization
    params.asl.do_normalization = true;
    params.asl.normvox = [2 2 2];

    % Smoothing
    params.asl.do_smoothing = true;
    params.asl.smoothfwhm = 6;

%% BE CAREFUL WITH CHANGING THE CODE BELOW THIS LINE !!
%---------------------------------------------------------------------------------------

restoredefaultpath

[params.my_spmbatch_path,~,~] = fileparts(mfilename('fullpath'));

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

my_spmbatch_start_aslpreprocessing(sublist,nsessions,datpath,params)

spm_figure('close',allchild(0));

cd(curdir)

if isfile(new_spm_read_vols_file), movefile(new_spm_read_vols_file,old_spm_read_vols_file); end

fprintf('\nDone\n')

end