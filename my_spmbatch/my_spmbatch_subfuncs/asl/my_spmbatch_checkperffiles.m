function ppparams = my_spmbatch_checkperffiles(ppparams,params)

namefilters(1).name = ppparams.substring;
namefilters(1).required = true;

namefilters(2).name = ppparams.sesstring;
namefilters(2).required = false;

namefilters(3).name = ['run-' num2str(ppparams.run)];
if params.func.mruns, namefilters(3).required = true; else namefilters(3).required = false; end

namefilters(4).name = ['task-' ppparams.task];
namefilters(4).required = true;

%% asl data

namefilters(5).name = '_asl';
namefilters(5).required = true;

aslniilist = my_spmbatch_dirfilelist(ppparams.subperfdir,'nii',namefilters,false);

if isempty(aslniilist)
    fprintf(['No asl nifti file found for ' ppparams.substring ' ' ppparams.sesstring ' task-' ppparams.task '\n'])
    fprintf('\nPP_Error\n');
    return
end

tmp = find(contains({aslniilist.name},'_echo-1'));
if ~isempty(tmp), aslniilist = aslniilist(tmp); 
else 
    tmp = find(contains({aslniilist.name(1:2)},'cd'));
    if ~isempty(tmp), aslniilist = aslniilist(tmp); end
end

prefixlist = split({aslniilist.name},'sub-');
if numel(aslniilist)==1, prefixlist=prefixlist(1); else prefixlist = prefixlist(:,:,1); end

for ip=1:numel(prefixlist)
    pref = prefixlist{ip};
    if contains(params.asl.splitaslbold,'meica'), test = strcmp(pref(1),'d'); end
    if contains(params.asl.splitaslbold,'dune'), test = strcmp(pref(1:2),'cd'); end

    if test
        ppparams.perf(1).aslprefix = pref; 
        fsplit = split(aslniilist(ip).name,ppparams.perf(1).aslprefix);
        ppparams.perf(1).aslfile = fsplit{2};
    end
end

%% label data

namefilters(5).name = '_label';
namefilters(5).required = true;

namefilters(6).name = '_echo-1';
namefilters(6).required = true;

labelniilist = my_spmbatch_dirfilelist(ppparams.subperfdir,'nii',namefilters,false);

if isempty(labelniilist)
    fprintf(['No label nifti file found for ' ppparams.substring ' ' ppparams.sesstring ' task-' ppparams.task '\n'])
    fprintf('\nPP_Error\n');
    return
end

prefixlist = split({labelniilist.name},'sub-');
if numel(labelniilist)==1, prefixlist=prefixlist(1); else prefixlist = prefixlist(:,:,1); end

tmp = find(contains(prefixlist,'f'));

if ~isempty(tmp)
    ppparams.perf(1).labelprefix = prefixlist{tmp};

    ffile = labelniilist(tmp).name;
    fsplit = split(ffile,ppparams.perf(1).labelprefix);
    ppparams.perf(1).labelfile = fsplit{end};
end

%% m0scan data

namefilters(5).name = '_m0scan';
namefilters(5).required = true;

namefilters(6).name = '_echo-1';
namefilters(6).required = true;

m0scanniilist = my_spmbatch_dirfilelist(ppparams.subperfdir,'nii',namefilters,false);

if isempty(m0scanniilist)
    fprintf(['No m0scan nifti file found for ' ppparams.substring ' ' ppparams.sesstring ' task-' ppparams.task '\n'])
    fprintf('\nPP_Error\n');
    return
end

prefixlist = split({m0scanniilist.name},'sub-');
if numel(m0scanniilist)==1, prefixlist=prefixlist(1); else prefixlist = prefixlist(:,:,1); end

fpresplit = split(ppparams.perf(1).aslprefix,'f');
studyprefix = fpresplit{end};

perfcheck = true;
while perfcheck
    tmp = find(strcmp(prefixlist,studyprefix));
    if ~isempty(tmp)
        ppparams.perf(1).m0scanprefix = studyprefix; 
        perfcheck = false;
    else 
        studyprefix = studyprefix(2:end); 
        if length(studyprefix) == 0 
            perfcheck = false; 
        end
    end
end

if ~isempty(tmp)
    ffile = m0scanniilist(tmp).name;
    fsplit = split(ffile,ppparams.perf(1).m0scanprefix);
    ppparams.perf(1).m0scanfile = fsplit{2};
end

%% Segmentation maps

if contains(params.asl.GMWM,'M0asl') 
    tmp=find(contains(prefixlist,'c1'));
    if ~isempty(tmp), ppparams.perf(1).c1m0scanfile = m0scanniilist(tmp).name; end
    
    tmp=find(contains(prefixlist,'c2'));
    if ~isempty(tmp), ppparams.perf(1).c2m0scanfile = m0scanniilist(tmp).name; end
    
    tmp=find(contains(prefixlist,'c3'));
    if ~isempty(tmp), ppparams.perf(1).c3m0scanfile = m0scanniilist(tmp).name; end
else
    ananamefilters(1).name = ppparams.substring;
    ananamefilters(1).required = true;
    
    ananamefilters(2).name = ppparams.sesstring;
    ananamefilters(2).required = false;
    
    ananamefilters(3).name = ['_T1w'];
    ananamefilters(3).required = true;

    anatniilist = my_spmbatch_dirfilelist(ppparams.subperfdir,'nii',ananamefilters,false);
    
    if ~isempty(anatniilist)
        prefixlist = split({anatniilist.name},'sub-');
        if numel(anatniilist)==1, prefixlist=prefixlist{1}; else prefixlist = prefixlist(:,:,1); end
        
        tmp=find(contains(prefixlist,'rp1e'));
        if ~isempty(tmp), ppparams.perf(1).c1m0scanfile = anatniilist(tmp).name; end
        
        tmp=find(contains(prefixlist,'rp2e'));
        if ~isempty(tmp), ppparams.perf(1).c2m0scanfile = anatniilist(tmp).name; end
        
        tmp=find(contains(prefixlist,'rp3e'));
        if ~isempty(tmp), ppparams.perf(1).c3m0scanfile = anatniilist(tmp).name; end
    end
end

%% deltam data

namefilters(5).name = '_deltam';
namefilters(5).required = true;

deltamniilist = my_spmbatch_dirfilelist(ppparams.subperfdir,'nii',namefilters,false);

if ~isempty(deltamniilist)
    for i=1:numel(deltamniilist)
        delete(fullfile(ppparams.subperfdir,deltamniilist(i).name))
    end
end

%% CBF data

namefilters(5).name = '_cbf';
namefilters(5).required = true;

cbfniilist = my_spmbatch_dirfilelist(ppparams.subperfdir,'nii',namefilters,false);

if ~isempty(cbfniilist)
    for i=1:numel(cbfniilist)
        delete(fullfile(ppparams.subperfdir,cbfniilist(i).name))
    end
end
