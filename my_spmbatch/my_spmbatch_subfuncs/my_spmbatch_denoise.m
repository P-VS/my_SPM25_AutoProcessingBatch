function [delfiles,keepfiles] = my_spmbatch_denoise(sub,ses,run,task,datpath,params)

delfiles = {};
keepfiles = {};

if ~params.func.meepi, params.func.echoes = [1]; end

ppparams.sub = sub;
ppparams.ses = ses;
ppparams.run = run;
ppparams.task = task;
ppparams.use_parallel = params.use_parallel;
ppparams.save_intermediate_results = params.save_intermediate_results;

%% Search for the data folders

ppparams.substring = ['sub-' num2str(sub,['%0' num2str(params.sub_digits) 'd'])];

ppparams.sesstring = ['ses-' num2str(ses,'%02d')];
if ~isfolder(fullfile(datpath,ppparams.substring,ppparams.sesstring)), ppparams.sesstring = ['ses-' num2str(ses,'%03d')]; end

ppparams.subpath = fullfile(datpath,ppparams.substring,ppparams.sesstring);

if ~isfolder(ppparams.subpath), ppparams.subpath = fullfile(datpath,ppparams.substring); end

if ~isfolder(ppparams.subpath)
    fprintf(['No data folder for subject ' num2str(sub) ' session ' num2str(ses)])
    fprintf('\nPP_Error\n');
    return
end

ppparams.subfuncdir = fullfile(ppparams.subpath,params.func_save_folder);

if ~isfolder(ppparams.subfuncdir)
    fprintf(['No preprocessed func data folder found for subject ' num2str(sub) ' session ' num2str(ses)])
    fprintf('\nPP_Error\n');
    return
end

ppparams = my_spmbatch_checkfuncfiles(ppparams,params,false);

if ~isfield(ppparams,'func') || isempty(ppparams.func(1).prefix), return; end
if contains(ppparams.func(1).prefix,'c')
    ppparams.func = ppparams.func(1);
    ppparams.echoes = 1;
    ppparams.meepi = false;
end

%% Denoising the ME/SE fMRI data
if ~contains(ppparams.func(1).prefix,'d')
    [ppparams,delfiles,keepfiles] = my_spmbatch_fmridenoising(ppparams,params,delfiles,keepfiles);

    keepfiles{numel(keepfiles)+1} = {fullfile(ppparams.subfuncdir,['d' ppparams.func(ie).prefix ppparams.func(ie).funcfile])};
end

%% Smooth denoised func
if ~contains(ppparams.func(1).prefix,'s')
    fprintf('Do smoothing \n')

    [ppparams,delfiles,keepfiles] = my_spmbatch_dosmoothfunc(ppparams,params,1,delfiles,keepfiles);
end