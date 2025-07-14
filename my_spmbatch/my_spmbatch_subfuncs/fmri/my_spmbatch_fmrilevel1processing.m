function matlabbatch = my_spmbatch_fmrilevel1processing(sub,ses,run,task,datpath,params)

matlabbatch = {};

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

ppparams.preprocfmridir = fullfile(ppparams.subpath,params.preprocfmridir);

if ~isfolder(ppparams.subfuncdir)
    fprintf(['No preprocessed func data folder found for subject ' num2str(sub) ' session ' num2str(ses)])
    fprintf('\nPP_Error\n');
    return
end

for ir=1:numel(params.iruns)
    %% Search for the data files
    
    namefilters(1).name = ppparams.substring;
    namefilters(1).required = true;
    
    namefilters(2).name = ppparams.sesstring;
    namefilters(2).required = false;
    
    switch params.func.use_runs
        case 'separately'
            namefilters(3).name = ['run-' num2str(run)];
        case 'together'
            namefilters(3).name = ['run-' num2str(params.iruns(ir))];
    end
    if params.func.mruns, namefilters(3).required = true; else namefilters(3).required = false; end
    
    namefilters(4).name = ['task-' task];
    namefilters(4).required = true;
    
    %fMRI data
    
    if contains(params.modality,'fmri'), namefilters(5).name = '_bold'; end
    if contains(params.modality,'fasl'), namefilters(5).name = '_cbf'; end
    namefilters(5).required = true;

    namefilters(6).name = params.fmri_prefix;
    namefilters(6).required = true;

    funcniilist = my_spmbatch_dirfilelist(ppparams.preprocfmridir,'nii',namefilters,false);
    
    if isempty(funcniilist)
        fprintf(['No nii files found for ' ppparams.substring ' ' ppparams.sesstring ' task-' params.task '\n'])
        fprintf('\nPP_Error\n');
        return
    end

    for ie=1:numel(params.func.echoes)
        if params.func.meepi && params.use_echoes_as_sessions %Filter list based on echo number
            tmp = find(or(contains({funcniilist.name},['_echo-' num2str(params.func.echoes(ie))]),contains({funcniilist.name},['_e' num2str(params.func.echoes(ie))])));
            if isempty(tmp), edirniilist = funcniilist; else edirniilist = funcniilist(tmp); end
        else
            edirniilist = funcniilist;
        end
    
        prefixlist = split({edirniilist.name},'sub-');
        if numel(edirniilist)==1, prefixlist=prefixlist{1}; else prefixlist = prefixlist(:,:,1); end
    
        tmp = find(strcmp(prefixlist,params.fmri_prefix));
        if ~isempty(tmp), ppparams.frun(ir).func(ie).funcfile = edirniilist(tmp).name; end
    
        if ~isfield(ppparams.frun(ir).func(ie),'funcfile')
            fprintf(['no preprocessed fmri data found for run ' num2str(ir) ' for echo ' num2str(params.func.echoes(ie)) '\n'])
            fprintf('\nPP_Error\n');
            return
        end
    
        Vfunc = spm_vol(fullfile(ppparams.preprocfmridir,ppparams.frun(ir).func(ie).funcfile));
    
        for i=1:numel(Vfunc)
            ppfmridat{ir}.sess{ie}.func{i,1} = [Vfunc(i).fname ',' num2str(i)];
        end
    end
    
    %confound file
    
    if params.add_regressors
        cnamefilters = namefilters;

        namefilters(5).name = '_bold';
        
        cnamefilters(6).name = params.confounds_prefix;
        cnamefilters(6).required = true;
    
        funcconlist = my_spmbatch_dirfilelist(ppparams.preprocfmridir,'txt',cnamefilters,false);
        
        if isempty(funcconlist)
            cnamefilters(5).name = '_asl';
            cnamefilters(5).required = true;
    
            funcconlist = my_spmbatch_dirfilelist(ppparams.preprocfmridir,'txt',cnamefilters,false);
        end

        if isempty(funcconlist)
            fprintf(['No confound file found for ' ppparams.substring ' ' ppparams.sesstring ' task-' task '\n'])
            fprintf('\nPP_Error\n');
            return
        end
        
        prefixlist = split({funcconlist.name},'sub-');
        if numel(funcconlist)==1, prefixlist=prefixlist{1}; else prefixlist = prefixlist(:,:,1); end
        
        tmp = find(strcmp(prefixlist,params.confounds_prefix));
        if ~isempty(tmp), ppparams.frun(ir).confoundsfile = fullfile(funcconlist(tmp).folder,funcconlist(tmp).name); else ppparams.frun(ir).confoundsfile = ''; end        
    else
        ppparams.frun(ir).confoundsfile = '';
    end

    %events.tsv file
    
    enamefilters(1:4) = namefilters(1:4);
    
    enamefilters(5).name = '_events';
    enamefilters(5).required = true;
    
    functsvlist = my_spmbatch_dirfilelist(ppparams.subfuncdir,'tsv',enamefilters,false);
    
    if isempty(functsvlist)
        fprintf(['No events.tsv files found for ' ppparams.substring ' ' ppparams.sesstring ' task-' task '\n'])
        fprintf('\nPP_Error\n');
        return
    end
    
    ppparams.frun(ir).functsvfile = fullfile(functsvlist(1).folder,functsvlist(1).name);
    
    %json file

    jnamefilters(1:4) = namefilters(1:4);
    
    jnamefilters(5).name = '_bold';
    jnamefilters(5).required = true;
    
    funcjsonlist = my_spmbatch_dirfilelist(ppparams.subfuncdir,'json',jnamefilters,false);
        
    if isempty(funcjsonlist)
        jnamefilters(5).name = '_asl';
        jnamefilters(5).required = true;
        
        funcjsonlist = my_spmbatch_dirfilelist(ppparams.subfuncdir,'json',jnamefilters,false);
    end

    if params.func.meepi
        jstmp = find(or(contains({funcjsonlist.name},'echo-1'),contains({funcjsonlist.name},'_e1')));
        funcjsonlist = funcjsonlist(jstmp);
    end
    
    if isempty(funcjsonlist)
        fprintf(['No json files found for ' ppparams.substring ' ' ppparams.sesstring ' task-' task '\n'])
        fprintf('\nPP_Error\n');
        return
    end
    
    ppparams.frun(ir).funcjsonfile = fullfile(funcjsonlist(1).folder,funcjsonlist(1).name);
end

%% fMRI model specification
if params.func.mruns && contains(params.func.use_runs,'separately')
    resultmap = fullfile(ppparams.subpath,['SPMMAT-' task '_' params.analysisname '_run-' num2str(run)]);
else
    resultmap = fullfile(ppparams.subpath,['SPMMAT-' task '_' params.analysisname]);
end

if exist(resultmap,'dir'); rmdir(resultmap,'s'); end
mkdir(resultmap)

%% Model Specification

jsondat = fileread(ppparams.frun(1).funcjsonfile);
jsondat = jsondecode(jsondat);

if contains(params.modality,'fasl') && params.reduced_temporal_resolution, tr=params.newTR; else tr = jsondat.RepetitionTime; end

if isfield(jsondat,'SliceTiming')
    SliceTiming = jsondat.SliceTiming;
    nsl= ceil(numel(SliceTiming)/numel(find(SliceTiming==SliceTiming(1))));
else
    nsl = 1;
end

matlabbatch{1}.spm.stats.fmri_spec.dir = {resultmap};
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = tr;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = nsl;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 1;

for ir=1:numel(params.iruns)
    % correct events file for dummy scans if needed
    dummys = floor(params.dummytime/tr);
    try
        edat{ir} = tdfread(ppparams.frun(ir).functsvfile,'\t');
    catch
        T = readtable(ppparams.frun(ir).functsvfile,'FileType','text');
        edat{ir}.onset = T.onset;
        edat{ir}.duration = T.duration;
        edat{ir}.trial_type = T.trial_type;
    end
    edat{ir}.onset = edat{ir}.onset-dummys*tr;
    
    for it=1:numel(edat{ir}.trial_type(:,1))
        ntrial_type(it,1) = convertCharsToStrings(edat{ir}.trial_type(it,:));
    end
    
    edat{ir}.trial_type = ntrial_type;
    
    [~,edatorder] = sort(edat{ir}.trial_type);
    edat{ir}.onset = edat{ir}.onset(edatorder);
    edat{ir}.duration = edat{ir}.duration(edatorder);
    edat{ir}.trial_type = edat{ir}.trial_type(edatorder,:);
    
    numc=0;
    for trial=1:numel(edat{ir}.onset)
        if isstring(edat{ir}.trial_type(trial,:))
            trial_type = edat{ir}.trial_type(trial,:);
        else
            trial_type = num2str(edat{ir}.trial_type(trial,:));
        end
    
        if numc>0
            numc = numel(edat{ir}.conditions)+1;
            for nc=1:numel(edat{ir}.conditions)
                if strcmp(edat{ir}.conditions{nc}.name,strtrim(trial_type)); numc=nc; end
            end
            if numc<numel(edat{ir}.conditions)+1
                edat{ir}.conditions{numc}.onsets = [edat{ir}.conditions{numc}.onsets edat{ir}.onset(trial)];
                edat{ir}.conditions{numc}.durations = [edat{ir}.conditions{numc}.durations edat{ir}.duration(trial)];
            else
                edat{ir}.conditions{numc}.name = strtrim(trial_type);
                edat{ir}.conditions{numc}.onsets = [edat{ir}.onset(trial)];
                edat{ir}.conditions{numc}.durations = [edat{ir}.duration(trial)];
            end
        else
            edat{ir}.conditions{1}.name = strtrim(trial_type);
            edat{ir}.conditions{1}.onsets = edat{ir}.onset(trial);
            edat{ir}.conditions{1}.durations = edat{ir}.duration(trial);
            numc=1;
        end
    end
    
    for ne=1:numel(params.func.echoes)
        nsess = (ir-1)*numel(params.func.echoes)+ne;
    
        matlabbatch{1}.spm.stats.fmri_spec.sess(nsess).scans = ppfmridat{ir}.sess{ne}.func(:,1);
    
        for nc=1:numel(edat{1}.conditions)
            matlabbatch{1}.spm.stats.fmri_spec.sess(nsess).cond(nc).name = char(edat{ir}.conditions{nc}.name);
            matlabbatch{1}.spm.stats.fmri_spec.sess(nsess).cond(nc).onset = edat{ir}.conditions{nc}.onsets;
            matlabbatch{1}.spm.stats.fmri_spec.sess(nsess).cond(nc).duration = edat{ir}.conditions{nc}.durations;
            matlabbatch{1}.spm.stats.fmri_spec.sess(nsess).cond(nc).tmod = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess(nsess).cond(nc).pmod = struct('name', {}, 'param', {}, 'poly', {});
            if params.reduced_temporal_resolution, matlabbatch{1}.spm.stats.fmri_spec.sess(nsess).cond(nc).orth = 0; 
            else matlabbatch{1}.spm.stats.fmri_spec.sess(nsess).cond(nc).orth = 1; end
        end
    
        matlabbatch{1}.spm.stats.fmri_spec.sess(nsess).multi = {''};

        matlabbatch{1}.spm.stats.fmri_spec.sess(nsess).regress = struct('name', {}, 'val', {});
    
        matlabbatch{1}.spm.stats.fmri_spec.sess(nsess).multi_reg = {ppparams.frun(ir).confoundsfile};

        if ~contains(ppparams.frun(1).func(1).funcfile,'f'), matlabbatch{1}.spm.stats.fmri_spec.sess(nsess).hpf = params.hpf;
        else
            Vfunc = spm_vol(fullfile(ppparams.preprocfmridir,ppparams.frun(1).func(1).funcfile));
            matlabbatch{1}.spm.stats.fmri_spec.sess(nsess).hpf = tr * (numel(Vfunc)-1);
        end
    end
end

matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
if ~params.add_derivatives; matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0]; 
else matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [1 1]; end
matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
matlabbatch{1}.spm.stats.fmri_spec.global = 'None';

Vfunc = spm_vol(fullfile(ppparams.preprocfmridir,ppparams.frun(1).func(1).funcfile));
nvols = min([numel(Vfunc),50]);
fdata = spm_read_vols(Vfunc(1:nvols));
mask = my_spmbatch_mask(fdata);

Vmask = Vfunc(1);
rmfield(Vmask,'pinfo');
Vmask.fname = fullfile(ppparams.preprocfmridir,['mask_' ppparams.frun(1).func(1).funcfile]);
Vmask.descrip = 'my_spmbatch - mask';
Vmask.dt = [spm_type('float32'),spm_platform('bigend')];
Vmask.n = [1 1];
Vmask = spm_write_vol(Vmask,mask);

mask_file = Vmask.fname;

clear fdata mask

if contains(params.modality,'fmri')
    matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.0;
elseif contains(params.modality,'fasl')
    matlabbatch{1}.spm.stats.fmri_spec.mthresh = -1.0; 
end
matlabbatch{1}.spm.stats.fmri_spec.cvi = params.model_serial_correlations;
matlabbatch{1}.spm.stats.fmri_spec.mask = {Vmask.fname};

%% Optimize GLM with TEDM
if params.optimize_HRF
    spm_jobman('run', matlabbatch)

    clear matlabbatch

    SPM_file = fullfile(resultmap,'SPM.mat');
    SPM_file = my_spmmbatch_tedm(SPM_file,resultmap,mask_file);

    mbidx = 1;
    matlabbatch{mbidx}.spm.stats.fmri_est.spmmat(1) = {SPM_file};
else
    mbidx = 2;
    matlabbatch{mbidx}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{mbidx-1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
end

%% Model estimation

matlabbatch{mbidx}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{mbidx}.spm.stats.fmri_est.method.Classical = 1;
  
mbidx=mbidx+1;

%% Contrast Manager

matlabbatch{mbidx}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{mbidx-1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));

for ic=1:numel(params.contrast)
    contrastname='';
    weights = [];

    for ir=1:numel(params.iruns)
        if ~params.add_derivatives; ncondcol = 1; else ncondcol = 3; end

        subweights = zeros(1,ncondcol*numel(numel(edat{ir}.onset)));

        if params.add_regressors
            rpdat = load(ppparams.frun(ir).confoundsfile);
            subweights=[subweights zeros(1,numel(rpdat(1,:)))];
        end
        for icn=1:numel(params.contrast(ic).conditions)
            if params.contrast(ic).vector(icn)>0; contrastname = [contrastname ' + ' params.contrast(ic).conditions{icn}]; end
            if params.contrast(ic).vector(icn)<0; contrastname = [contrastname ' - ' params.contrast(ic).conditions{icn}]; end
        
            indx=0;
            for icn2=1:numel(edat{ir}.conditions)
                if strcmp(lower(params.contrast(ic).conditions{icn}),lower(edat{ir}.conditions{icn2}.name)); indx=(icn2-1)*ncondcol+1; end
            end
    
            if indx>0; subweights(indx:indx+(ncondcol-1))=params.contrast(ic).vector(icn); end
        end

        subweights = repmat(subweights,1,numel(params.func.echoes));

        weights=[weights subweights];
    end

    matlabbatch{mbidx}.spm.stats.con.consess{ic}.tcon.name = contrastname;
    matlabbatch{mbidx}.spm.stats.con.consess{ic}.tcon.weights = weights;
    
    matlabbatch{mbidx}.spm.stats.con.consess{ic}.tcon.sessrep = 'none';
end

matlabbatch{mbidx}.spm.stats.con.delete = 0;

mbidx=mbidx+1;

%% SPM Results
if params.save_spm_results
    matlabbatch{mbidx}.spm.stats.results.spmmat = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{mbidx-1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{mbidx}.spm.stats.results.conspec.titlestr = '';
    matlabbatch{mbidx}.spm.stats.results.conspec.contrasts = Inf;
    matlabbatch{mbidx}.spm.stats.results.conspec.threshdesc = params.threshold_correction;
    matlabbatch{mbidx}.spm.stats.results.conspec.thresh = params.pthreshold;
    matlabbatch{mbidx}.spm.stats.results.conspec.extent = params.kthreshold;
    matlabbatch{mbidx}.spm.stats.results.conspec.conjunction = 1;
    matlabbatch{mbidx}.spm.stats.results.conspec.mask.none = 1;
    matlabbatch{mbidx}.spm.stats.results.units = 1;

    oi = 1;
    if params.save_thresholded_map
        matlabbatch{mbidx}.spm.stats.results.export{oi}.tspm.basename = 'thres';
        oi = oi+1;
    end
    if params.save_binary_mask
        matlabbatch{mbidx}.spm.stats.results.export{oi}.binary.basename = 'bin';
        oi = oi+1;
    end
    if params.save_naray
        matlabbatch{mbidx}.spm.stats.results.export{oi}.nary.basename = 'n-aray';
        oi = oi+1;
    end
    if params.save_csv_file
        matlabbatch{mbidx}.spm.stats.results.export{oi}.csv = true;
        oi = oi+1;
    end
    if params.save_pdf_file
        matlabbatch{mbidx}.spm.stats.results.export{oi}.pdf = true;
        oi = oi+1;
    end
    if params.save_tiff_file
        matlabbatch{mbidx}.spm.stats.results.export{oi}.tif = true;
        oi = oi+1;
    end
end
