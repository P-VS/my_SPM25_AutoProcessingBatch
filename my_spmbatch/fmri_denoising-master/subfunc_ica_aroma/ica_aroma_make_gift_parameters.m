function icaparms_file = ica_aroma_make_gift_parameters(ica_source_file,t_r,n_sessions,mask,comp_est,ppparams)

icaparms_file = fullfile(ppparams.subfuncdir,'input_spatial_ica.m');

fid = fopen(icaparms_file,'w');

fprintf(fid,'%s','modalityType = ''fMRI'';');

%% Type of stability analysis
% Options are 1 and 2.
% 1 - Regular Group ICA
% 2 - Group ICA using icasso
% 3 - Group ICA using Minimum spanning tree (MST)
fprintf(fid,'\n%s','which_analysis = 1;');

%% ICASSO options.
% This variable will be used only when which_analysis variable is set to 2.
fprintf(fid,'\n%s','icasso_opts.sel_mode = ''bootstrap'';');  % Options are 'randinit', 'bootstrap' and 'both'
fprintf(fid,'\n%s','icasso_opts.num_ica_runs = 5;'); % Number of times ICA will be run
% Most stable run estimate is based on these settings. 
fprintf(fid,'\n%s','icasso_opts.min_cluster_size = 4;'); % Minimum cluster size
fprintf(fid,'\n%s','icasso_opts.max_cluster_size = 5;'); % Max cluster size. Max is the no. of components

%% Enter TR in seconds.
fprintf(fid,'\n%s',['TR = ' num2str(t_r) ';']);

%% Group ica type
% Options are spatial or temporal for fMRI modality. By default, spatial
% ica is run if not specified.
fprintf(fid,'\n%s','group_ica_type = ''spatial'';');

%% Parallel info
% enter mode serial or parallel. If parallel, enter number of
% sessions/workers to do job in parallel
v = ver;
if any(strcmp('Parallel Computing Toolbox', {v.Name})) %&& ~ppparams.use_parallel
    fprintf(fid,'\n%s','parallel_info.mode = ''parallel'';');
    fprintf(fid,'\n%s','parallel_info.num_workers = 4;');
else
    fprintf(fid,'\n%s','parallel_info.mode = ''serial'';');
end

%% Group PCA performance settings. Best setting for each option will be selected based on variable MAX_AVAILABLE_RAM in icatb_defaults.m. 
% If you have selected option 3 (user specified settings) you need to manually set the PCA options. See manual or other
% templates (icatb/icatb_batch_files/Input_data_subjects_1.m) for more information to set PCA options 
%
% Options are:
% 1 - Maximize Performance
% 2 - Less Memory Usage
% 3 - User Specified Settings
fprintf(fid,'\n%s','perfType = 1;');

%% Conserve disk space
% Conserve disk space. Options are:
% 0 - Write all analysis files including intermediate files (PCA, Backreconstruction, scaled component MAT files)
% 1 - Write only necessary files required to resume the analysis. The files written are as follows:
%   a. Data reduction files - Only eigen vectors and eigen values are written in the first data reduction step. PCA components are written at the last reduction stage.
%   b. Back-reconstruction files - Back-reconstruction files are not written to the disk. The information is computed while doing scaling components step
%   c. Scaling components files - Scaling components MAT files are not written when using GIFT or SBM.
%   d. Group stats files - Only mean of all data-sets is written.
% 2 - Write all files till the group stats. Cleanup intermediate files at the end of the group stats (PCA, Back-reconstruct, Scaled component MAT files in GIFT). Analysis cannot be
% resumed if there are any changes to the setup parameters. Utilities that work with PCA and Backreconstruction files like Remove components, Percent Variance, etc won't work with this option .
fprintf(fid,'\n%s','conserve_disk_space = 2;');

%% Design matrix selection
% Design matrix (SPM.mat) is used for sorting the components
% temporally (time courses) during display. Design matrix will not be used during the
% analysis stage except for SEMI-BLIND ICA.
% options are ('no', 'same_sub_same_sess', 'same_sub_diff_sess', 'diff_sub_diff_sess')
% 1. 'no' - means no design matrix.
% 2. 'same_sub_same_sess' - same design over subjects and sessions
% 3. 'same_sub_diff_sess' - same design matrix for subjects but different
% over sessions
% 4. 'diff_sub_diff_sess' - means one design matrix per subject.
fprintf(fid,'\n%s','keyword_designMatrix = ''no'';');

%% There are three ways to enter the subject data
% options are 1, 2, 3 or 4
fprintf(fid,'\n%s','dataSelectionMethod = 2;');

%% Method 2
% If you have different filePatterns and location for subjects not in one
% root folder then enter the data here.
% Number of subjects is determined getting the length of the selected subjects. specify the data set or data sets needed for 
% the analysis here.
fprintf(fid,'\n%s','selectedSubjects = {''s1''};');  % naming for subjects s1 refers to subject 1, Use cell array convention even in case of one subject one session

% Number of Sessions
fprintf(fid,'\n%s',['numOfSess = ' num2str(n_sessions) ';']);

% functional data folder, file pattern and file numbers to include
% You can provide the file numbers ([1:220]) to include as a vector. If you want to
% select all the files then leave empty.

[spth,snm,~] = fileparts(ica_source_file);

fprintf(fid,'\n%s',['s1_s1 = {']);
fprintf(fid,'%s',['''' spth ''',''' [snm '.nii'] '''']);
fprintf(fid,'%s','};');

%for i=1:numel(ppparams.echoes)
%    Vfunc = spm_vol(fullfile(ppparams.subfuncdir,[ppparams.func(i).prefix ppparams.func(i).funcfile]));
    
%    fprintf(fid,'\n%s',['s1_s' num2str(i) ' = {']);
%    fprintf(fid,'%s',['''' ppparams.subfuncdir ''',''' [ppparams.func(i).prefix ppparams.func(i).funcfile] '''']); % ',(1:' num2str(numel(Vfunc)) ')']); % subject 1 session 1
%    fprintf(fid,'%s','};');
%end

%% Enter directory to put results of analysis
outdir = fullfile(ppparams.subfuncdir, ['ica-dir_' ppparams.task]);
fprintf(fid,'\n%s',['outputDir = ''' outdir ''';']);

%% Enter Name (Prefix) Of Output Files
fprintf(fid,'\n%s','prefix = ''ica_aroma'';');

%% Enter location (full file path) of the image file to use as mask
% or use Default mask which is []
fprintf(fid,'\n%s',['maskFile = ''' mask ''';']);

%% Group PCA Type. Used for analysis on multiple subjects and sessions when 2 data reduction steps are used.
% Options are 'subject specific' and 'grand mean'. 
%   a. Subject specific - Individual PCA is done on each data-set before group
%   PCA is done.
%   b. Grand Mean - PCA is done on the mean over all data-sets. Each data-set is
%   projected on to the eigen space of the mean before doing group PCA.
%
% NOTE: Grand mean implemented is from FSL Melodic. Make sure that there are
% equal no. of timepoints between data-sets.
%
fprintf(fid,'\n%s','group_pca_type = ''Grand Mean'';');

%% Back reconstruction type. Options are 1 and 2
% 1 - Regular
% 2 - Spatial-temporal Regression 
% 3 - GICA3
% 4 - GICA
% 5 - GIG-ICA
fprintf(fid,'\n%s','backReconType = 4;');

%% Data Pre-processing options
% 1 - Remove mean per time point
% 2 - Remove mean per voxel
% 3 - Intensity normalization
% 4 - Variance normalization
fprintf(fid,'\n%s','preproc_type = 1;');

%% Maximum reduction steps you can select is 2
% You have the option to select one data-reduction or 2 data reduction
% steps when spatial ica is used. For temporal ica, only one data-reduction
% is done.
fprintf(fid,'\n%s','numReductionSteps = 1;');

%% Batch Estimation. If 1 is specified then estimation of 
% the components takes place and the corresponding PC numbers are associated
% Options are 1 or 0
fprintf(fid,'\n%s','doEstimation = 0;');

%% Number of pc to reduce each subject down to at each reduction step
% The number of independent components the will be extracted is the same as 
% the number of principal components after the final data reduction step.  
fprintf(fid,'\n%s',['numOfPC1 = ' num2str(comp_est) ';']);
fprintf(fid,'\n%s','numOfPC2 = 0;');

%% Scale the Results. Options are 0, 1, 2
% 0 - Don't scale
% 1 - Scale to Percent signal change
% 2 - Scale to Z scores
fprintf(fid,'\n%s','scaleType = 2;');

%% 'Which ICA Algorithm Do You Want To Use';
% see icatb_icaAlgorithm for details or type icatb_icaAlgorithm at the
% command prompt.
% Note: Use only one subject and one session for Semi-blind ICA. Also specify atmost two reference function names

% 1 means infomax, 2 means fastICA, etc.
fprintf(fid,'\n%s','algoType = 1;');

fclose(fid);