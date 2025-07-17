function ppparams = my_spmbatch_checkfuncfiles(ppparams,params,for_functional)

%% list of files

%functional nii files
namefilters(1).name = ppparams.substring;
namefilters(1).required = true;

namefilters(2).name = ppparams.sesstring;
namefilters(2).required = false;

namefilters(3).name = ['run-' num2str(ppparams.run)];
if params.func.mruns, namefilters(3).required = true; else namefilters(3).required = false; end

namefilters(4).name = ['task-' ppparams.task];
namefilters(4).required = true;

namefilters(5).name = '_bold';
namefilters(5).required = true;

funcniilist = my_spmbatch_dirfilelist(ppparams.subfuncdir,'nii',namefilters,false);

if params.func.isaslbold
    namefilters(5).name = '_aslbold';
    namefilters(5).required = true;
    
    aslniilist = my_spmbatch_dirfilelist(ppparams.subfuncdir,'nii',namefilters,false);

    if isempty(funcniilist), funcniilist=aslniilist; else funcniilist = [funcniilist;aslniilist]; end
end

if isempty(funcniilist)
    fprintf(['No nifti files found for ' ppparams.substring ' ' ppparams.sesstring ' task-' ppparams.task '\n'])
    fprintf('\nPP_Error\n');
    ppparams.error = true;
    return
end

tmp=find(~contains({funcniilist.name},'_e_'));
if ~isempty(tmp), funcniilist = funcniilist(tmp); end

%functional json files
funcjsonlist = my_spmbatch_dirfilelist(ppparams.subfuncdir,'json',namefilters,false);

if isempty(funcjsonlist)
    subjfuncdir = fullfile(ppparams.subpath,'func');

    funcjsonlist = my_spmbatch_dirfilelist(subjfuncdir,'json',namefilters,false);

    if isempty(funcjsonlist)
        fprintf(['No json files found for ' ppparams.substring ' ' ppparams.sesstring ' task-' ppparams.task '\n'])
        fprintf('\nPP_Error\n');
        ppparams.error = true;
        return
    end
end

%physiologic noise files (motion, derivatives, aCompCor)
funcrplist = my_spmbatch_dirfilelist(ppparams.subfuncdir,'txt',namefilters,false);

if isempty(funcrplist) && for_functional
    subtfuncdir = fullfile(ppparams.subpath,'func');

    funcrplist = my_spmbatch_dirfilelist(subtfuncdir,'txt',namefilters,false);
end

if ~isempty(funcrplist)
    rpprefixlist = split({funcrplist.name},'sub-');
    if numel(rpprefixlist)==1, rpprefixlist=rpprefixlist(1); else rpprefixlist = rpprefixlist(:,:,1); end
end

ppparams.echoes = params.func.echoes;
ppparams.meepi = params.func.meepi;

%% Check steps done on functional data
for ie=params.func.echoes
    if params.func.meepi %Filter list based on echo number
        tmp = find(or(contains({funcniilist.name},['_echo-' num2str(ie)]),contains({funcniilist.name},['_e' num2str(ie)])));
        if ~isempty(tmp), edirniilist = funcniilist(tmp); else edirniilist = []; end
    
        if ie==params.func.echoes(1)
            tmp = find(and(~contains({funcniilist.name},'_echo-'),~contains({funcniilist.name},'_e')));
            if ~isempty(tmp), edirniilist = [edirniilist;funcniilist(tmp)]; end
        end

        jstmp = find(or(contains({funcjsonlist.name},['echo-' num2str(ie)]),contains({funcjsonlist.name},['_e' num2str(ie)])));
        if ~isempty(jstmp), ppparams.func(ie).jsonfile = fullfile(funcjsonlist(jstmp(1)).folder,funcjsonlist(jstmp(1)).name); end
    else
        edirniilist = funcniilist;

        ppparams.func(ie).jsonfile = fullfile(funcjsonlist(1).folder,funcjsonlist(1).name);
    end

    if ~isempty(edirniilist)
        prefixlist = split({edirniilist.name},'sub-');
        if numel(edirniilist)==1, prefixlist=prefixlist(1); else prefixlist = prefixlist(:,:,1); end
    
        if ~for_functional
            %% For denoising after normalization
            studyprefix = 'sdw';
            ppparams.func(ie).prefix = '';

            perfcheck = true;
            while perfcheck
                tmp = find(and(contains(prefixlist,studyprefix),~contains(prefixlist,['s' studyprefix])));
                if ~isempty(tmp)
                    ppparams.func(ie).prefix = prefixlist{tmp}; 
                    perfcheck = false;
                else 
                    studyprefix = studyprefix(2:end); 
                    if length(studyprefix) == 0
                        perfcheck = false; 
                    end
                end
            end

            if ~isempty(tmp)
                ffile = edirniilist(tmp).name;
                fsplit = split(ffile,ppparams.func(ie).prefix);
                ppparams.func(ie).funcfile = fsplit{2};
            end
        else
            %% For preprocessing
            tmp = find(strlength(prefixlist)==0);
            if ~isempty(tmp), ppparams.func(ie).funcfile = edirniilist(tmp).name; end
        
            ppparams.func(ie).prefix = '';
            studyprefix = 'e';
        
            if params.func.do_realignment, studyprefix = ['r' studyprefix]; end
            if params.func.pepolar, studyprefix = ['u' studyprefix]; end
            if params.func.isaslbold, studyprefix = ['f' studyprefix]; end
            if params.func.do_ArtRepair, studyprefix = ['v' studyprefix]; end
            if params.func.do_slicetime, studyprefix = ['a' studyprefix]; end
            if params.func.isaslbold && params.denoise.do_DUNE, studyprefix = ['l' studyprefix]; end
            if ~params.func.isaslbold && params.func.denoise && params.denoise.do_bpfilter, studyprefix = ['f' studyprefix]; end
            if params.func.denoise && (params.denoise.do_noiseregression || params.denoise.do_ICA_AROMA), studyprefix = ['d' studyprefix]; end
            if params.func.meepi && params.func.do_echocombination, studyprefix = ['c' studyprefix]; end
            if params.func.do_normalization, studyprefix = ['w' studyprefix]; end
            if params.func.do_smoothing, studyprefix = ['s' studyprefix]; end
        
            perfcheck = true;
            while perfcheck
                tmp = find(strcmp(prefixlist,studyprefix));
                if ~isempty(tmp)
                    ppparams.func(ie).prefix = studyprefix; 
                    perfcheck = false;
                else 
                    studyprefix = studyprefix(2:end); 
                    if length(studyprefix) == 0; perfcheck = false; end
                end
            end
        
            if params.func.do_realignment && ~(contains(ppparams.func(ie).prefix,'r') || contains(ppparams.func(ie).prefix,'u'))
                if ~isfile(fullfile(ppparams.subfuncdir,['e' ppparams.func(params.func.echoes(1)).funcfile])), ppparams.func(params.func.echoes(1)).prefix = ''; end
            end
        end
    end
end

%% Check for extra files

%saved preprocessed results

namefilters(5).name = '_bold';
namefilters(5).required = true;

subsavedir = fullfile(ppparams.subpath,params.func_save_folder);
if isfolder(subsavedir) && ~for_functional, funcniilist = my_spmbatch_dirfilelist(subsavedir,'nii',namefilters,false); end

if params.func.meepi %Filter list based on echo number
    tmp = find(or(contains({funcniilist.name},'_echo-1'),contains({funcniilist.name},'_e1')));
    if ~isempty(tmp), edirniilist = funcniilist(tmp); end

    tmp = find(and(~contains({funcniilist.name},'_echo-'),~contains({funcniilist.name},'_e')));
    if ~isempty(tmp), edirniilist = [edirniilist;funcniilist(tmp)]; end
else
    edirniilist = funcniilist;
end

prefixlist = split({edirniilist.name},'sub-');
if numel(edirniilist)==1, prefixlist=prefixlist(1); else prefixlist = prefixlist(:,:,1); end

if params.func.do_realignment || params.do_denoising
    if ~isempty(funcrplist)
        %rp_... file
        tmp = find(and(contains(rpprefixlist,'rp_e'),and(~contains(rpprefixlist,'der_'),~contains(rpprefixlist,'acc_'))));
        if ~isempty(tmp), ppparams.rp_file = fullfile(funcrplist(tmp).folder,funcrplist(tmp).name); end
   
        %der_... file
        tmp = find(contains(rpprefixlist,'der_'));
        if ~isempty(tmp), ppparams.der_file = fullfile(funcrplist(tmp).folder,funcrplist(tmp).name); end

        %acc_... file
        tmp = find(contains(rpprefixlist,'acc_'));
        if ~isempty(tmp), ppparams.acc_file = fullfile(funcrplist(tmp).folder,funcrplist(tmp).name); end
    end

    denpreff = '';
    if contains(ppparams.func(params.func.echoes(1)).prefix,'e'), denpreff = ['e' denpreff]; end
    if contains(ppparams.func(params.func.echoes(1)).prefix,'r'), denpreff = ['r' denpreff]; end
    if contains(ppparams.func(params.func.echoes(1)).prefix,'u'), denpreff = ['u' denpreff]; end
    if params.func.isaslbold && contains(ppparams.func(params.func.echoes(1)).prefix,'f'), denpreff = ['f' denpreff]; end
    if contains(ppparams.func(params.func.echoes(1)).prefix,'v'), denpreff = ['v' denpreff]; end
    if contains(ppparams.func(params.func.echoes(1)).prefix,'a'), denpreff = ['a' denpreff]; end

    if params.do_denoising
        if contains(ppparams.func(params.func.echoes(1)).prefix,'c'), denpreff = ['c' denpreff]; end
        if contains(ppparams.func(params.func.echoes(1)).prefix,'w'), denpreff = ['w' denpreff]; end
    end

    %mask used for denoising if before normalization
    tmp = find(strcmp(prefixlist,['fmask_' denpreff]));
    if ~isempty(tmp), ppparams.fmask = fullfile(edirniilist(tmp).folder,edirniilist(tmp).name); end

    %segmentation used for denoising if before normalization
    tmp = find(strcmp(prefixlist,['c1' denpreff]));
    if ~isempty(tmp), ppparams.fc1im = fullfile(edirniilist(tmp).folder,edirniilist(tmp).name); end

    tmp = find(strcmp(prefixlist,['c2' denpreff]));
    if ~isempty(tmp), ppparams.fc2im = fullfile(edirniilist(tmp).folder,edirniilist(tmp).name); end

    tmp = find(strcmp(prefixlist,['c3' denpreff]));
    if ~isempty(tmp), ppparams.fc3im = fullfile(edirniilist(tmp).folder,edirniilist(tmp).name); end

    %acc_... file
    if ~isempty(funcrplist)
        tmp = find(strcmp(rpprefixlist,['acc_' denpreff]));
        if ~isempty(tmp), ppparams.acc_file = fullfile(funcrplist(tmp).folder,funcrplist(tmp).name); end
    end

    if ~isempty(funcrplist)
        %ica_... file
        tmp = find(strcmp(rpprefixlist,['nboldica_' denpreff]));
        if ~isempty(tmp), ppparams.nboldica_file = fullfile(funcrplist(tmp).folder,funcrplist(tmp).name); end

        tmp = find(strcmp(rpprefixlist,['naslica_' denpreff]));
        if ~isempty(tmp), ppparams.naslica_file = fullfile(funcrplist(tmp).folder,funcrplist(tmp).name); end
    end
end

%% fmap data for geometric correction
if for_functional
    if params.func.pepolar
        
        fmfilters(1).name = ppparams.substring;
        fmfilters(1).required = true;
        
        fmfilters(2).name = ppparams.sesstring;
        fmfilters(2).required = false;
        
        fmfilters(3).name = ['run-' num2str(ppparams.run)];
        fmfilters(3).required = false;
     
        fmfilters(4).name = 'dir-';
        fmfilters(4).required = true;

        fmfilters(4).name = '_epi';
        fmfilters(4).required = true;
    
        fmapniilist = my_spmbatch_dirfilelist(ppparams.subfmapdir,'nii',fmfilters,true);
        fmapjsonlist = my_spmbatch_dirfilelist(ppparams.subfmapdir,'json',fmfilters,true);
    
        if isempty(fmapniilist) || isempty(fmapjsonlist)
            params.func.pepolar = false;
    
            fprintf(['No fmap files found for ' ppparams.substring ' ' ppparams.sesstring ' run ' num2str(ppparams.run) ' task-' ppparams.task '\n'])
        end

        tmp=find(~contains({fmapniilist.name},'_e_'));
        if ~isempty(tmp), fmapniilist = fmapniilist(tmp); end
    
        %% check for fmap data in case of pepolar
        for ie=params.func.echoes
            if params.func.meepi %Filter list based on echo number
                tmp = find(or(contains({fmapniilist.name},['echo-' num2str(ie)]),contains({fmapniilist.name},['_e' num2str(ie)])));
                if isempty(tmp)
                    fprintf(['no fmap data found for echo ' num2str(ie) '\n'])
                    fprintf('\nPP_Error\n');
                    ppparams.error = true;
                    return
                end
            
                efmapniilist = fmapniilist(tmp);
        
                jstmp = find(or(contains({fmapjsonlist.name},['echo-' num2str(ie)]),contains({fmapjsonlist.name},['_e' num2str(ie)])));
                if isempty(jstmp)
                    fprintf(['no fmap json file found for echo ' num2str(ie) '\n'])
                    fprintf('\nPP_Error\n');
                    ppparams.error = true;
                    return
                end
        
                ppparams.func(ie).fmapjsonfile = fullfile(fmapjsonlist(jstmp(1)).folder,fmapjsonlist(jstmp(1)).name);
    
            else
                efmapniilist = fmapniilist;
        
                ppparams.func(ie).fmapjsonfile = fullfile(fmapjsonlist(1).folder,fmapjsonlist(1).name);
            end
        
            fmprefixlist = split({efmapniilist.name},'sub-');
            fmprefixlist = fmprefixlist(:,:,1);
                
            tmp = find(strlength(fmprefixlist)==0);
            if ~isempty(tmp)
                ppparams.func(ie).fmapfile = efmapniilist(tmp).name; 
            else
                fprintf(['no fmap data found for echo ' num2str(ie) '\n'])
                fprintf('\nPP_Error\n');
                params.func.pepolar = false;
            end
    
            ppparams.func(ie).fmap_prefix = '';
        end
    end
end