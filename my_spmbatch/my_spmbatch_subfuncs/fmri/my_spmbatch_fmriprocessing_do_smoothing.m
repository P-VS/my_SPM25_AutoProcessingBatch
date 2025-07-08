function [params] = my_spmbatch_fmriprocessing_do_smoothing(sub,ses,run,task,datpath,params)

delfiles = {};
keepfiles = {};

params.save_folder = params.preprocfmridir;
params.func.smoothfwhm = 6; %(default=6)

if ~params.func.meepi, params.func.echoes = [1]; end

ppparams.sub = sub;
ppparams.ses = ses;
ppparams.run = run;
ppparams.task = task;
ppparams.use_parallel = params.use_parallel;
ppparams.save_intermediate_results = false;

%% Search for the data folders

ppparams.substring = ['sub-' num2str(sub,'%02d')];
if ~isfolder(fullfile(datpath,ppparams.substring)), ppparams.substring = ['sub-' num2str(sub,'%03d')]; end

ppparams.sesstring = ['ses-' num2str(ses,'%02d')];
if ~isfolder(fullfile(datpath,ppparams.substring,ppparams.sesstring)), ppparams.sesstring = ['ses-' num2str(ses,'%03d')]; end

ppparams.subpath = fullfile(datpath,ppparams.substring,ppparams.sesstring);

if ~isfolder(ppparams.subpath), ppparams.subpath = fullfile(datpath,ppparams.substring); end

if ~isfolder(ppparams.subpath)
    fprintf(['No data folder for subject ' num2str(sub) ' session ' num2str(ses)])
    fprintf('\nPP_Error\n');
    return
end

ppparams.subfuncdir = fullfile(ppparams.subpath,params.save_folder);

if ~isfolder(ppparams.subfuncdir)
    fprintf(['No preprocessed func data folder found for subject ' num2str(sub) ' session ' num2str(ses)])
    fprintf('\nPP_Error\n');
    return
end

%% Search for the data files

namefilters(1).name = ppparams.substring;
namefilters(1).required = true;

namefilters(2).name = ppparams.sesstring;
namefilters(2).required = false;

namefilters(3).name = ['run-' num2str(run)];
if params.func.mruns, namefilters(3).required = true; else namefilters(3).required = false; end

namefilters(4).name = ['task-' task];
namefilters(4).required = true;

namefilters(5).name = '_bold';
namefilters(5).required = true;

namefilters(6).name = params.fmri_prefix;
namefilters(6).required = true;

funcniilist = my_spmbatch_dirfilelist(ppparams.subfuncdir,'nii',namefilters,false);

if isempty(funcniilist)
    fprintf(['No nifti files found for ' substring ' ' sesstring ' task-' task '\n'])
    fprintf('\nPP_Error\n');
    return
end

%% do smoothing
for ie=params.func.echoes
    if params.func.meepi %Filter list based on echo number
        tmp = find(or(contains({funcniilist.name},['_echo-' num2str(ie)]),contains({funcniilist.name},['_e' num2str(ie)])));
        if isempty(tmp), funcfile = funcniilist(1).name; else funcfile = funcniilist(tmp).name; end
    else
        funcfile = funcniilist(1).name;
    end
    ppparams.func(ie).prefix = params.fmri_prefix;

    ffile = funcfile;
    fsplit = split(ffile,ppparams.func(ie).prefix);
    ppparams.func(ie).funcfile = fsplit{2};

    [ppparams,~,~] = my_spmbatch_dosmoothfunc(ppparams,params,ie,delfiles,keepfiles);
end

params.fmri_prefix = ppparams.func(1).prefix;