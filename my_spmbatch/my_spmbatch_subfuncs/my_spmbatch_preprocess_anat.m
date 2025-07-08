function [delfiles,keepfiles] = my_spmbatch_preprocess_anat(sub,ses,datpath,params)

delfiles = {};
keepfiles = {};

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

ppparams.subanatdir = fullfile(ppparams.subpath,'anat');

if ~isfolder(ppparams.subanatdir)
    fprintf(['No anat data folder found for subject ' num2str(sub) ' session ' num2str(ses)])
    fprintf('\nPP_Error\n');
    return
end

ppparams.preproc_anat = fullfile(ppparams.subpath,'preproc_anat');
if ~exist(ppparams.preproc_anat,"dir"), mkdir(ppparams.preproc_anat); end

%% Search for the data files

namefilters(1).name = ppparams.substring;
namefilters(1).required = true;

namefilters(2).name = ppparams.sesstring;
namefilters(2).required = false;

namefilters(3).name = ['_T1w'];
namefilters(3).required = true;

% Unprocessed data

anatniilist = my_spmbatch_dirfilelist(ppparams.subanatdir,'nii',namefilters,false);

if isempty(anatniilist)
    fprintf(['No nifti files found for ' ppparams.substring ' ' ppparams.sesstring '\n'])
    fprintf('\nPP_Error\n');
    return
end

tmp = find(contains({anatniilist.name},'_Crop_1'));
if ~isempty(tmp), anatniilist = anatniilist(tmp); end

prefixlist = split({anatniilist.name},'sub-');
if numel(anatniilist)==1, prefixlist=prefixlist{1}; else prefixlist = prefixlist(:,:,1); end

tmp = find(strlength(prefixlist)==0);
if ~isempty(tmp), ppparams.subanat = anatniilist(tmp).name; end

% Normalized data

ppanatniilist = my_spmbatch_dirfilelist(ppparams.preproc_anat,'nii',namefilters,false);

if ~isempty(ppanatniilist)
    tmp = find(contains({ppanatniilist.name},'_Crop_1'));
    if ~isempty(tmp), ppanatniilist = ppanatniilist(tmp); end
    
    prefixlist = split({ppanatniilist.name},'sub-');
    prefixlist = prefixlist(:,:,1);
end
    
anatprefix = 'e';

if ~isempty(ppanatniilist) && params.anat.do_normalization
    tmp = find(strcmp(prefixlist,['w' anatprefix]));

    if isempty(tmp), tmp = find(strcmp(prefixlist,['wm' anatprefix])); end

    if ~isempty(tmp), ppparams.wsubanat = ppanatniilist(tmp).name; end
end

% Segmented data

if ~isempty(ppanatniilist) && params.anat.do_segmentation
    tmp = find(strcmp(prefixlist,['wp1' anatprefix]));
    if isempty(tmp), tmp = find(strcmp(prefixlist,['wc1' anatprefix])); end

    if ~isempty(tmp), ppparams.anat.wc1im = panatniilist(tmp).name; end

    tmp = find(strcmp(prefixlist,['wp2' anatprefix]));
    if isempty(tmp), tmp = find(strcmp(prefixlist,['wc2' anatprefix])); end

    if ~isempty(tmp), ppparams.anat.wc2im = panatniilist(tmp).name; end

    tmp = find(strcmp(prefixlist,['wp3' anatprefix]));
    if isempty(tmp), tmp = find(strcmp(prefixlist,['wc3' anatprefix])); end

    if ~isempty(tmp), ppparams.anat.wc3im = panatniilist(tmp).name; end 
end

%% Do reoriention

nm = split(ppparams.subanat,'.nii');
transfile = fullfile(ppparams.subanatdir,[nm{1} '_reorient.mat']);
if isfile(transfile)
    load(transfile,'M')
    transM = M;
else        
    transM = my_spmbatch_vol_set_com(fullfile(ppparams.subanatdir,ppparams.subanat));
    transM(1:3,4) = -transM(1:3,4);
end

Vanat = spm_vol(fullfile(ppparams.subanatdir,ppparams.subanat));
MM = Vanat.private.mat0;

Vanat = my_reset_orientation(Vanat,transM * MM);

anatdat = spm_read_vols(Vanat);

Vanat.fname = fullfile(ppparams.subanatdir,['e' ppparams.subanat]);
Vanat.descrip = 'reoriented';
Vanat = spm_create_vol(Vanat);
Vanat = spm_write_vol(Vanat,anatdat);

auto_acpc_reorient(Vanat.fname,'T1');

ppparams.subanat = ['e' ppparams.subanat];
delfiles{numel(delfiles)+1} = {fullfile(ppparams.subanatdir,ppparams.subanat)};

%% Do segmentation

if params.anat.do_segmentation && ~isfield(ppparams,'wc1im') && ~isfield(ppparams,'wc2im') && ~isfield(ppparams,'wc3im')
    fprintf('Do segmentation \n')

    % Normalization
    params.vbm.do_normalization = params.anat.do_normalization;
    params.vbm.normvox = params.anat.normvox;

    % Segmentation
    params.vbm.do_segmentation = true;
    if isfield (params.anat,'roi_atlas'), params.vbm.do_roi_atlas = params.anat.roi_atlas; else params.vbm.do_roi_atlas = false; end
    if isfield (params.anat,'do_surface'), params.vbm.do_surface = params.anat.do_surface; else params.vbm.do_surface = false; end

    [delfiles,keepfiles] = my_spmbatch_cat12vbm(ppparams,params,delfiles,keepfiles);

    ppparams.anat.wc1im = ['mwp1' ppparams.subanat];
    ppparams.anat.wc2im = ['mwp2' ppparams.subanat];
    ppparams.anat.wc3im = ['mwp3' ppparams.subanat];
    ppparams.wsubanat = ['wm1' ppparams.subanat];

    if ~params.anat.do_normalization
        ppparams.anat.c1im = ['p1' ppparams.subanat];
        ppparams.anat.c2im = ['p2' ppparams.subanat];
        ppparams.anat.c3im = ['p3' ppparams.subanat];
    end
end

%% Do Normalization

if params.anat.do_normalization && ~isfield(ppparams,'wsubanat')
    fprintf('Do normalization \n')

    anatnormestwrite.subj.vol = {fullfile(ppparams.subanatdir,ppparams.subanat)};
    anatnormestwrite.subj.resample = {fullfile(ppparams.subanatdir,ppparams.subanat)};
    anatnormestwrite.eoptions.biasreg = 0.0001;
    anatnormestwrite.eoptions.biasfwhm = 60;
    anatnormestwrite.eoptions.tpm = {fullfile(spm('Dir'),'tpm','TPM.nii')};
    anatnormestwrite.eoptions.affreg = 'mni';
    anatnormestwrite.eoptions.reg = [0 0.001 0.5 0.05 0.2];
    anatnormestwrite.eoptions.fwhm = 0;
    anatnormestwrite.eoptions.samp = 3;
    anatnormestwrite.woptions.bb = [-78 -112 -70;78 76 85];
    anatnormestwrite.woptions.vox = params.anat.normvox;
    anatnormestwrite.woptions.interp = 4;
    anatnormestwrite.woptions.prefix = 'w';

    spm_run_norm(anatnormestwrite);

    ppparams.deffile = fullfile(ppparams.subanatdir,['y_' ppparams.subanat]);
    delfiles{numel(delfiles)+1} = {ppparams.deffile};
    keepfiles{numel(keepfiles)+1} = {fullfile(ppparams.subanatdir,['w' ppparams.subanat])};

    ppparams.wsubanat = ['w' ppparams.subanat];
end