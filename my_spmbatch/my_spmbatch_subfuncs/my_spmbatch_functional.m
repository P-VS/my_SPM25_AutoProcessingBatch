function [delfiles,keepfiles] = my_spmbatch_functional(sub,ses,run,task,datpath,params)

delfiles = {};
keepfiles = {};

if ~params.func.meepi, params.func.echoes = [1]; end

ppparams.sub = sub;
ppparams.ses = ses;
ppparams.run = run;
ppparams.task = task;
ppparams.echoes = params.func.echoes;
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

ppparams.subfuncdir = fullfile(ppparams.subpath,'func');

if ~isfolder(ppparams.subfuncdir)
    fprintf(['No func data folder found for subject ' num2str(sub) ' session ' num2str(ses)])
    fprintf('\nPP_Error\n');
    return
end

if params.func.pepolar
    ppparams.subfmapdir = fullfile(ppparams.subpath,'fmap');

    if ~isfolder(ppparams.subfmapdir)
        fprintf(['No fmap data folder found for subject ' num2str(sub) ' session ' num2str(ses)])
        fprintf('\nPP_Error\n');

        ppparams.subfmapdir = ppparams.subfuncdir;
    end
end

%% Search for the data files
ppparams = my_spmbatch_checkfuncfiles(ppparams,params,true);

%% Do the preprocessing
[~,delfiles,keepfiles] = my_spmbatch_preprocfunc(ppparams,params,delfiles,keepfiles);