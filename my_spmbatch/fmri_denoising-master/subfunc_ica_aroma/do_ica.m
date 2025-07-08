function do_ica(ica_source_file,mask, t_r, ppparams)
% DO_ICA Creates the parameter file necessary to run ICA with the GIFT toolbox
% and performs ICA

%if numel(ppparams.echoes)>1 && ~contains(ppparams.func(1).prefix,'c')
%    n_sessions = numel(ppparams.echoes);
%else
%    n_sessions = 1;
%end

%% Estimate number of components (MDL - FWHM)
dim_est_opts.method = 2; %MDL
dim_est_opts.fwhm = [5,5,5]; %FWHM
comp_est = 0;

[comp_est, ~, ~, ~] = icatb_estimate_dimension(ica_source_file, mask, 'double', dim_est_opts);  
if comp_est>100, comp_est=100; end

%% Adapt template parameter to specific subject data

icatb_defaults;

% Load template parameters .mat file
icaparms_file = ica_aroma_make_gift_parameters(ica_source_file,t_r,1,mask,comp_est,ppparams);

%% Set up ICA
icaparam_file = icatb_setup_analysis(icaparms_file);

load(icaparam_file);

%% Run Analysis (All steps)
sesInfo = icatb_runAnalysis(sesInfo, [2:6]);

end


