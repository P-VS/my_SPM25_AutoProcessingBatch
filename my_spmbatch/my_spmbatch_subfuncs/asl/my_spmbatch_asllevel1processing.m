function matlabbatch = my_spmbatch_asllevel1processing(sub,ses,run,task,datpath,params)

params.add_parametricModulation = false; %use the weights in events.tsv for parametric modultion
params.add_regressors = false; %if data not denoised set true otherwhise false 
params.add_derivatives = false; %add temmperal and dispertion derivatives to the GLM (default=false)
params.optimize_HRF = false;

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
    
    %fASL data
    
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

    if params.func.meepi && params.use_echoes_as_sessions %Filter list based on echo number
        tmp = find(or(contains({funcniilist.name},['_echo-' num2str(params.func.echoes(1))]),contains({funcniilist.name},['_e' num2str(params.func.echoes(1))])));
        if isempty(tmp), edirniilist = funcniilist; else edirniilist = funcniilist(tmp); end
    else
        edirniilist = funcniilist;
    end

    prefixlist = split({edirniilist.name},'sub-');
    if numel(edirniilist)==1, prefixlist=prefixlist{1}; else prefixlist = prefixlist(:,:,1); end

    tmp = find(strcmp(prefixlist,params.fmri_prefix));
    if ~isempty(tmp), ppparams.frun(ir).func(1).funcfile = edirniilist(tmp).name; end

    if ~isfield(ppparams.frun(ir).func(1),'funcfile')
        fprintf(['no preprocessed fmri data found for run ' num2str(ir) ' for echo ' num2str(params.func.echoes(1)) '\n'])
        fprintf('\nPP_Error\n');
        return
    end

    Vfunc = spm_vol(fullfile(ppparams.preprocfmridir,ppparams.frun(ir).func(1).funcfile));

    for i=1:numel(Vfunc)
        ppfmridat{ir}.sess{1}.func{i,1} = [Vfunc(i).fname ',' num2str(i)];
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
    
    jnamefilters(5).name = '_asl';
    jnamefilters(5).required = true;
    
    funcjsonlist = my_spmbatch_dirfilelist(ppparams.subfuncdir,'json',jnamefilters,false);
        
    jstmp = find(or(contains({funcjsonlist.name},'echo-1'),contains({funcjsonlist.name},'_e1')));
    funcjsonlist = funcjsonlist(jstmp);
    
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

matlabbatch{1}.spm.stats.factorial_design.dir = {resultmap};

%% Model Specification

jsondat = fileread(ppparams.frun(1).funcjsonfile);
jsondat = jsondecode(jsondat);

tr = jsondat.RepetitionTime;

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
    edat{ir}.onset = edat{ir}.onset-(params.asl.LabelingDuration+params.asl.PostLabelDelay);
    
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

    if numel(edat{1}.conditions)==1
        edat{ir}.conditions{2}.name = 'rest';
        if ~edat{ir}.conditions{1}.onsets(1)==0
            edat{ir}.conditions{2}.onsets = 0;
            edat{ir}.conditions{2}.durations = edat{ir}.conditions{1}.onsets(1);
        end

        for iblock=1:numel(edat{ir}.conditions{1}.onsets)
            if edat{ir}.conditions{1}.onsets(iblock)+edat{ir}.conditions{1}.durations(iblock)<(numel(ppfmridat{ir}.sess{1}.func)*tr)
                nonset = edat{ir}.conditions{1}.onsets(iblock)+edat{ir}.conditions{1}.durations(iblock);
                if iblock<numel(edat{ir}.conditions{1}.onsets)
                    nduration = edat{ir}.conditions{1}.onsets(iblock+1)-nonset;
                else
                    nduration = numel(ppfmridat{ir}.sess{1}.func)*tr-nonset;
                end
                if iblock==1 && edat{ir}.conditions{1}.onsets(1)==0
                    edat{ir}.conditions{2}.onsets = nonset;
                    edat{ir}.conditions{2}.durations = nduration;
                else
                    edat{ir}.conditions{2}.onsets = [edat{ir}.conditions{2}.onsets nonset];
                    edat{ir}.conditions{2}.durations = [edat{ir}.conditions{2}.durations nduration];
                end
            end
        end

        for ic=1:numel(params.contrast)
            params.contrast(ic).conditions = {params.contrast(ic).conditions{1} 'rest'};
            if params.contrast(ic).vector(1)>0, params.contrast(ic).vector = [params.contrast(1).vector -1];
            else params.contrast(ic).vector = [params.contrast(ic).vector 1]; end
        end
    end

    for nc=1:numel(edat{1}.conditions)
        edat{ir}.conditions{nc}.startblock=floor(1+edat{ir}.conditions{nc}.onsets/tr);
        edat{ir}.conditions{nc}.endblock=floor((edat{ir}.conditions{nc}.onsets+edat{ir}.conditions{nc}.durations)/tr); 

        for iblock=1:numel(edat{ir}.conditions{nc}.onsets)
            if iblock==1
                edat{ir}.conditions{nc}.scans = [edat{ir}.conditions{nc}.startblock(iblock):edat{ir}.conditions{nc}.endblock(iblock)];
            else
                edat{ir}.conditions{nc}.scans = [edat{ir}.conditions{nc}.scans,edat{ir}.conditions{nc}.startblock(iblock):edat{ir}.conditions{nc}.endblock(iblock)];
            end
        end
        matlabbatch{1}.spm.stats.factorial_design.des.anova.icell(nc).scans = ppfmridat{ir}.sess{1}.func(edat{ir}.conditions{nc}.scans,1);
    end
end

matlabbatch{1}.spm.stats.factorial_design.des.anova.dept = 1;
matlabbatch{1}.spm.stats.factorial_design.des.anova.variance = 0;
matlabbatch{1}.spm.stats.factorial_design.des.anova.gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.anova.ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});

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

matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 0;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {Vmask.fname};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

%% Model estimation

matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

%% Contrast Manager

matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));

for ic=1:numel(params.contrast)
    contrastname='';
    weights = [];

    for ir=1:numel(params.iruns)
        if ~params.add_derivatives; ncondcol = 1; else ncondcol = 3; end
        subweights = zeros(1,ncondcol*numel(edat{ir}.conditions));
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
    
            if indx>0; subweights(indx)=params.contrast(ic).vector(icn); end
        end

        subweights = repmat(subweights,1,numel(params.func.echoes));

        weights=[weights subweights];
    end

    matlabbatch{3}.spm.stats.con.consess{ic}.tcon.name = contrastname;
    matlabbatch{3}.spm.stats.con.consess{ic}.tcon.weights = weights;
    
    matlabbatch{3}.spm.stats.con.consess{ic}.tcon.sessrep = 'none';
end

matlabbatch{3}.spm.stats.con.delete = 0;

%% SPM Results
if params.save_spm_results
    matlabbatch{4}.spm.stats.results.spmmat = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{4}.spm.stats.results.conspec.titlestr = '';
    matlabbatch{4}.spm.stats.results.conspec.contrasts = Inf;
    matlabbatch{4}.spm.stats.results.conspec.threshdesc = params.threshold_correction;
    matlabbatch{4}.spm.stats.results.conspec.thresh = params.pthreshold;
    matlabbatch{4}.spm.stats.results.conspec.extent = params.kthreshold;
    matlabbatch{4}.spm.stats.results.conspec.conjunction = 1;
    matlabbatch{4}.spm.stats.results.conspec.mask.none = 1;
    matlabbatch{4}.spm.stats.results.units = 1;

    oi = 1;
    if params.save_thresholded_map
        matlabbatch{4}.spm.stats.results.export{oi}.tspm.basename = 'thres';
        oi = oi+1;
    end
    if params.save_binary_mask
        matlabbatch{4}.spm.stats.results.export{oi}.binary.basename = 'bin';
        oi = oi+1;
    end
    if params.save_naray
        matlabbatch{4}.spm.stats.results.export{oi}.nary.basename = 'n-aray';
        oi = oi+1;
    end
    if params.save_csv_file
        matlabbatch{4}.spm.stats.results.export{oi}.csv = true;
        oi = oi+1;
    end
    if params.save_pdf_file
        matlabbatch{4}.spm.stats.results.export{oi}.pdf = true;
        oi = oi+1;
    end
    if params.save_tiff_file
        matlabbatch{4}.spm.stats.results.export{oi}.tif = true;
        oi = oi+1;
    end
end