function [delfiles,keepfiles] = my_spmbatch_aslpreprocessed(sub,ses,run,datpath,params)

delfiles = {};
keepfiles = {};

ppparams.sub = sub;
ppparams.ses = ses;
ppparams.run = run;
ppparams.reorient = params.reorient;

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

ppparams.subasldir = fullfile(ppparams.subpath,'perf');

if ~isfolder(ppparams.subasldir)
    fprintf(['No perf data folder found for subject ' num2str(sub) ' session ' num2str(ses)])
    fprintf('\nPP_Error\n');
    return
end

%% Search for the data files

namefilters(1).name = ppparams.substring;
namefilters(1).required = true;

namefilters(2).name = ppparams.sesstring;
namefilters(2).required = false;

namefilters(3).name = ['run-' num2str(ppparams.run)];
if params.asl.mruns, namefilters(3).required = true; else namefilters(3).required = false; end

namefilters(4).name = ['_deltam'];
namefilters(4).required = true;

%deltam

deltamniilist = my_spmbatch_dirfilelist(ppparams.subasldir,'nii',namefilters,false);

if isempty(deltamniilist)
    fprintf(['No deltam nifti files found for ' ppparams.substring ' ' ppparams.sesstring '\n'])
    fprintf('\nPP_Error\n');
    return
end

prefixlist = split({deltamniilist.name},'sub-');
if numel(deltamniilist)==1, prefixlist=prefixlist{1}; else prefixlist = prefixlist(:,:,1); end

tmp = find(strlength(prefixlist)==0);
if ~isempty(tmp), ppparams.deltam = fullfile(deltamniilist(tmp).folder,deltamniilist(tmp).name); end

tmp = find(strcmp(prefixlist,'e'));
if ~isempty(tmp), ppparams.edeltam = fullfile(deltamniilist(tmp).folder,deltamniilist(tmp).name); end

tmp = find(strcmp(prefixlist,'we'));
if ~isempty(tmp), ppparams.wdeltam = fullfile(deltamniilist(tmp).folder,deltamniilist(tmp).name); end

% deltam json

deltamjsonlist = my_spmbatch_dirfilelist(ppparams.subasldir,'json',namefilters,false);

if isempty(deltamjsonlist)
    fprintf(['No json files found for ' ppparams.substring ' ' ppparams.sesstring '\n'])
    fprintf('\nPP_Error\n');
    return
end

ppparams.deltamjson = fullfile(deltamjsonlist(1).folder,deltamjsonlist(1).name);

% m0scan

namefilters(4).name = ['_m0scan'];
namefilters(4).required = true;

m0scanniilist = my_spmbatch_dirfilelist(ppparams.subasldir,'nii',namefilters,false);

if isempty(m0scanniilist)
    fprintf(['No m0scan nifti files found for ' ppparams.substring ' ' ppparams.sesstring '\n'])
    fprintf('\nPP_Error\n');
    return
end

prefixlist = split({m0scanniilist.name},'sub-');
if numel(m0scanniilist)==1, prefixlist=prefixlist{1}; else prefixlist = prefixlist(:,:,1); end

tmp = find(strlength(prefixlist)==0);
if ~isempty(tmp), ppparams.m0scan = fullfile(m0scanniilist(tmp).folder,m0scanniilist(tmp).name); end

tmp = find(strcmp(prefixlist,'e'));
if ~isempty(tmp), ppparams.em0scan = fullfile(m0scanniilist(tmp).folder,m0scanniilist(tmp).name); end

tmp = find(strcmp(prefixlist,'re'));
if ~isempty(tmp), ppparams.rm0scan = fullfile(m0scanniilist(tmp).folder,m0scanniilist(tmp).name); end

tmp = find(strcmp(prefixlist,'tre'));
if ~isempty(tmp), ppparams.tm0scan = fullfile(m0scanniilist(tmp).folder,m0scanniilist(tmp).name); end

tmp = find(strcmp(prefixlist,'wre'));
if ~isempty(tmp), ppparams.wm0scan = fullfile(m0scanniilist(tmp).folder,m0scanniilist(tmp).name); end

% m0scan json

m0scanjsonlist = my_spmbatch_dirfilelist(ppparams.subasldir,'json',namefilters,false);

if isempty(m0scanjsonlist)
    fprintf(['No json files found for ' ppparams.substring ' ' ppparams.sesstring '\n'])
    fprintf('\nPP_Error\n');
    return
end

ppparams.m0scanjson = fullfile(m0scanjsonlist(1).folder,m0scanjsonlist(1).name);

% cbf map

namefilters(4).name = ['_cbf'];
namefilters(4).required = true;

cbfniilist = my_spmbatch_dirfilelist(ppparams.subasldir,'nii',namefilters,false);

if ~isempty(cbfniilist)
    prefixlist = split({cbfniilist.name},'sub-');
    if numel(edirniilist)==1, prefixlist=prefixlist{1}; else prefixlist = prefixlist(:,:,1); end

    tmp = find(strlength(prefixlist)==0);
    if ~isempty(tmp), ppparams.cbfmap = fullfile(cbfniilist(tmp).folder,cbfniilist(tmp).name); end
    
    tmp = find(strcmp(prefixlist,'e'));
    if ~isempty(tmp), ppparams.ecbfmap = fullfile(cbfniilist(tmp).folder,cbfniilist(tmp).name); end
    
    tmp = find(strcmp(prefixlist,'qe'));
    if ~isempty(tmp), ppparams.qcbfmap = fullfile(cbfniilist(tmp).folder,cbfniilist(tmp).name); end
    
    tmp = find(strcmp(prefixlist,'wqe'));
    if ~isempty(tmp), ppparams.wcbfmap = fullfile(cbfniilist(tmp).folder,cbfniilist(tmp).name); end
    
    tmp = find(strcmp(prefixlist,'swqe'));
    if ~isempty(tmp), ppparams.scbfmap = fullfile(cbfniilist(tmp).folder,cbfniilist(tmp).name); end
end

if ~isfield(ppparams,'cbfmap'), params.asl.do_cbfmapping = true; end

[delfiles,keepfiles] = my_spmbatch_preprocasl_processed(ppparams,params,delfiles,keepfiles);