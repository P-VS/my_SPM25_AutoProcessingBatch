function [ppparams,delfiles,keepfiles] = my_spmbatch_asl_segmentation(ppparams,params,delfiles,keepfiles)

if contains(params.asl.GMWM,'M0asl')
    %% Do segmentation on M0 image
    preproc.channel.vols = {fullfile(ppparams.subperfdir,[ppparams.perf(1).m0scanprefix ppparams.perf(1).m0scanfile]);};
    preproc.channel.biasreg = 0.001;
    preproc.channel.biasfwhm = 60;
    preproc.channel.write = [0 0];
    preproc.tissue(1).tpm = {fullfile(spm('Dir'),'tpm','TPM.nii,1')};
    preproc.tissue(1).ngaus = 1;
    preproc.tissue(1).native = [1 0];
    preproc.tissue(1).warped = [0 0];
    preproc.tissue(2).tpm = {fullfile(spm('Dir'),'tpm','TPM.nii,2')};
    preproc.tissue(2).ngaus = 1;
    preproc.tissue(2).native = [1 0];
    preproc.tissue(2).warped = [0 0];
    preproc.tissue(3).tpm = {fullfile(spm('Dir'),'tpm','TPM.nii,3')};
    preproc.tissue(3).ngaus = 2;
    preproc.tissue(3).native = [1 0];
    preproc.tissue(3).warped = [0 0];
    preproc.tissue(4).tpm = {fullfile(spm('Dir'),'tpm','TPM.nii,4')};
    preproc.tissue(4).ngaus = 3;
    preproc.tissue(4).native = [0 0];
    preproc.tissue(4).warped = [0 0];
    preproc.tissue(5).tpm = {fullfile(spm('Dir'),'tpm','TPM.nii,5')};
    preproc.tissue(5).ngaus = 4;
    preproc.tissue(5).native = [0 0];
    preproc.tissue(5).warped = [0 0];
    preproc.tissue(6).tpm = {fullfile(spm('Dir'),'tpm','TPM.nii,6')};
    preproc.tissue(6).ngaus = 2;
    preproc.tissue(6).native = [0 0];
    preproc.tissue(6).warped = [0 0];
    preproc.warp.mrf = 1;
    preproc.warp.cleanup = 1;
    preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
    preproc.warp.affreg = 'mni';
    preproc.warp.fwhm = 0;
    preproc.warp.samp = 3;
    preproc.warp.write = [0 1];
    preproc.warp.vox = NaN;
    preproc.warp.bb = [NaN NaN NaN;NaN NaN NaN];
    
    spm_preproc_run(preproc);
    
    ppparams.perf(1).c1m0scanfile = ['c1' ppparams.perf(1).m0scanprefix ppparams.perf(1).m0scanfile];
    ppparams.perf(1).c2m0scanfile = ['c2' ppparams.perf(1).m0scanprefix ppparams.perf(1).m0scanfile];
    ppparams.perf(1).c3m0scanfile = ['c3' ppparams.perf(1).m0scanprefix ppparams.perf(1).m0scanfile];
    
    fname = split([ppparams.perf(1).m0scanprefix ppparams.perf(1).m0scanfile],'.nii');
     
    delfiles{numel(delfiles)+1} = {fullfile(ppparams.subperfdir,[fname{1} '._seg8.mat'])};
    delfiles{numel(delfiles)+1} = {fullfile(ppparams.subperfdir,['y_' fname{1} '.nii'])};  
    delfiles{numel(delfiles)+1} = {fullfile(ppparams.subperfdir,ppparams.perf(1).c1m0scanfile)};
    delfiles{numel(delfiles)+1} = {fullfile(ppparams.subperfdir,ppparams.perf(1).c2m0scanfile)};
    delfiles{numel(delfiles)+1} = {fullfile(ppparams.subperfdir,ppparams.perf(1).c3m0scanfile)};

else
    %% Search for the anatomical data files

    ppparams.subanatdir = fullfile(ppparams.subpath,'anat');
    
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

    params.vbm.do_normalization = false;
    params.vbm.normvox = params.anat.normvox;

    %% Do Segmentation
    params.vbm.do_segmentation = true;
    params.vbm.do_roi_atlas = false;
    params.vbm.do_surface = false;

    [delfiles,~] = my_spmbatch_cat12vbm(ppparams,params,delfiles,keepfiles);

    ppparams.anat.c1im = ['p1' ppparams.subanat];
    ppparams.anat.c2im = ['p2' ppparams.subanat];
    ppparams.anat.c3im = ['p3' ppparams.subanat];
    ppparams.subanat = ['p0' ppparams.subanat];

    delfiles{numel(delfiles)+1} = {fullfile(ppparams.subanatdir,ppparams.subanat)};
    delfiles{numel(delfiles)+1} = {fullfile(ppparams.subanatdir,ppparams.anat.c1im)};
    delfiles{numel(delfiles)+1} = {fullfile(ppparams.subanatdir,ppparams.anat.c2im)};
    delfiles{numel(delfiles)+1} = {fullfile(ppparams.subanatdir,ppparams.anat.c3im)};

    %% Do Coregistration to M0 image
     
    estwrite.ref(1) = {fullfile(ppparams.subperfdir,[ppparams.perf(1).m0scanprefix ppparams.perf(1).m0scanfile])};
    estwrite.source(1) = {fullfile(ppparams.subanatdir,ppparams.subanat)};
    estwrite.other = {fullfile(ppparams.subanatdir,ppparams.anat.c1im), ...
                      fullfile(ppparams.subanatdir,ppparams.anat.c2im), ...
                      fullfile(ppparams.subanatdir,ppparams.anat.c3im)};
    estwrite.eoptions.cost_fun = 'nmi';
    estwrite.eoptions.sep = [4 2];
    estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    estwrite.eoptions.fwhm = [7 7];
    estwrite.roptions.interp = 4;
    estwrite.roptions.wrap = [0 0 0];
    estwrite.roptions.mask = 0;
    estwrite.roptions.prefix = 'r';
    
    out_coreg = spm_run_coreg(estwrite);
    
    ppparams.anat.c1im = ['r' ppparams.anat.c1im];
    ppparams.anat.c2im = ['r' ppparams.anat.c2im];
    ppparams.anat.c3im = ['r' ppparams.anat.c3im];

    %% Copy to perfussion map

    movefile(fullfile(ppparams.subanatdir,ppparams.anat.c1im),fullfile(ppparams.subperfdir,ppparams.anat.c1im))
    movefile(fullfile(ppparams.subanatdir,ppparams.anat.c2im),fullfile(ppparams.subperfdir,ppparams.anat.c2im))
    movefile(fullfile(ppparams.subanatdir,ppparams.anat.c3im),fullfile(ppparams.subperfdir,ppparams.anat.c3im))

    ppparams.perf(1).c1m0scanfile = ppparams.anat.c1im;
    ppparams.perf(1).c2m0scanfile = ppparams.anat.c2im;
    ppparams.perf(1).c3m0scanfile = ppparams.anat.c3im;
end
