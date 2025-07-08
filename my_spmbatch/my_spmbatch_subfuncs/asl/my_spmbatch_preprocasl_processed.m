function [delfiles,keepfiles] = my_spmbatch_preprocasl_processed(ppparams,params,delfiles,keepfiles)

%% Load and reorient data

if params.reorient && ~isfield(ppparams,'em0scan') && ~isfield(ppparams,'edeltam')
    fprintf('Start reorient \n')

    %m0scan

    [pth nm ext] = fileparts(ppparams.m0scan);
    transfile = fullfile(ppparams.subasldir,[nm '_reorient.mat']);
    if isfile(transfile)
        load(transfile,'M')
        transM = M;
    else        
        transM = my_spmbatch_vol_set_com(ppparams.m0scan);
        transM(1:3,4) = -transM(1:3,4);
    end

    Vm0 = spm_vol(ppparams.m0scan);
    MM = Vm0.private.mat0;

    Vm0 = my_reset_orientation(Vm0,transM*MM);

    m0dat = spm_read_vols(Vm0);

    m0mask = my_spmbatch_mask(m0dat);
    m0dat(m0mask<0.5) = 0;

    Vm0.fname = spm_file(ppparams.m0scan, 'prefix','e');
    Vm0.descrip = 'reoriented';
    Vm0 = spm_create_vol(Vm0);
    Vm0 = spm_write_vol(Vm0,m0dat);

    ppparams.em0scan = spm_file(ppparams.m0scan, 'prefix','e');
    delfiles{numel(delfiles)+1} = {ppparams.em0scan};

    auto_acpc_reorient(ppparams.em0scan,'PD');

    Vm0 = spm_vol(ppparams.em0scan);
    MM = Vm0.mat;

    %deltam

    Vmd = spm_vol(ppparams.deltam);

    Vmd = my_reset_orientation(Vmd,MM);

    mddat = spm_read_vols(Vmd);

    Vmd.fname = spm_file(ppparams.deltam, 'prefix','e');
    Vmd.descrip = 'reoriented';
    Vmd = spm_create_vol(Vmd);
    Vmd = spm_write_vol(Vmd,mddat);

    ppparams.edeltam = spm_file(ppparams.deltam, 'prefix','e');
    delfiles{numel(delfiles)+1} = {ppparams.edeltam};

    %cbfmap

    if isfield(ppparams,'cbfmap')
        Vcbf = spm_vol(ppparams.cbfmap);
    
        Vcbf = my_reset_orientation(Vcbf,MM);
    
        cbfdat = spm_read_vols(Vcbf);
    
        Vcbf.fname = spm_file(ppparams.cbfmap, 'prefix','e');
        Vcbf.descrip = 'reoriented';
        Vcbf = spm_create_vol(Vcbf);
        Vcbf = spm_write_vol(Vcbf,cbfdat);
    
        ppparams.ecbfmap = spm_file(ppparams.cbfmap, 'prefix','e');
        delfiles{numel(delfiles)+1} = {ppparams.ecbfmap};
    end
elseif ~isfield(ppparams,'em0scan') && ~isfield(ppparams,'edeltam')
    ppparams.em0scan = spm_file(ppparams.m0scan, 'prefix','e');

    copyfile(ppparams.m0scan,ppparams.em0scan);

    ppparams.edeltam = spm_file(ppparams.deltam, 'prefix','e');

    copyfile(ppparams.deltam,ppparams.edeltam);

    if isfield(ppparams,'cbfmap') && ~isfield(ppparams,'ecbfmap')
        ppparams.ecbfmap = spm_file(ppparams.cbfmap, 'prefix','e');

        copyfile(ppparams.cbfmap,ppparams.ecbfmap);
    end
end

%% Corregistration M0map to deltam
if ~isfield(ppparams,'rm0scan')
    fprintf('Start corregistration \n')
    
    estwrite.ref(1) = {ppparams.edeltam};
    estwrite.source(1) = {ppparams.em0scan};
    estwrite.other = {ppparams.em0scan};
    estwrite.eoptions.cost_fun = 'nmi';
    estwrite.eoptions.sep = [4 2];
    estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    estwrite.eoptions.fwhm = [7 7];
    estwrite.roptions.interp = 4;
    estwrite.roptions.wrap = [0 0 0];
    estwrite.roptions.mask = 0;
    estwrite.roptions.prefix = 'r';
    
    out_coreg = spm_run_coreg(estwrite);
    
    ppparams.rm0scan = spm_file(ppparams.em0scan, 'prefix','r');
    delfiles{numel(delfiles)+1} = {ppparams.rm0scan};
end

%% T1 correction of M0 and Quantification of CBF
if params.asl.do_cbfmapping || ~isfield(ppparams,'ecbfmap')
    if ~isfield(ppparams,'tm0scan')
        fprintf('Start M0 preprocessing \n')
    
        [ppparams,delfiles,keepfiles] = my_spmbatch_asl_M0correction(ppparams,params,delfiles,keepfiles);
    end

    if ~isfield(ppparams,'qcbfmap')
        fprintf('Start CBF quantification \n')
    
        [ppparams,delfiles,keepfiles] = my_spmbatch_asl_cbfquantification(ppparams,delfiles,keepfiles);
    end
elseif isfield(ppparams,'ecbfmap')
    ppparams.qcbfmap = ppparams.ecbfmap;
end

%% Normalization of the ASL scan
if params.asl.do_normalization
    if ~isfield(ppparams,'wm0scan') || ~isfield(ppparams,'wdeltam') || ~isfield(ppparams,'wcbfmap')
        fprintf('Start normalization \n')
    
        resamplescans = {ppparams.rm0scan,ppparams.edeltam,ppparams.qcbfmap};
    
        aslnormestwrite.subj.vol = {ppparams.rm0scan};
        aslnormestwrite.subj.resample = resamplescans;
        aslnormestwrite.eoptions.biasreg = 0.0001;
        aslnormestwrite.eoptions.biasfwhm = 60;
        aslnormestwrite.eoptions.tpm = {fullfile(spm('Dir'),'tpm','TPM.nii')};
        aslnormestwrite.eoptions.affreg = 'mni';
        aslnormestwrite.eoptions.reg = [0 0.001 0.5 0.05 0.2];
        aslnormestwrite.eoptions.fwhm = 0;
        aslnormestwrite.eoptions.samp = 3;
        aslnormestwrite.woptions.bb = [-78 -112 -70;78 76 85];
        aslnormestwrite.woptions.vox = params.asl.normvox;
        aslnormestwrite.woptions.interp = 1;
        aslnormestwrite.woptions.prefix = 'w';
    
        spm_run_norm(aslnormestwrite);
    
        ppparams.deffile = spm_file(ppparams.rm0scan, 'prefix','y_');
        delfiles{numel(delfiles)+1} = {ppparams.deffile};
    
        keepfiles{numel(keepfiles)+1} = {spm_file(ppparams.rm0scan, 'prefix','w')};
        keepfiles{numel(keepfiles)+1} = {spm_file(ppparams.edeltam, 'prefix','w')};
    
        if isfield(ppparams,'qcbfmap')
            keepfiles{numel(keepfiles)+1} = {spm_file(ppparams.qcbfmap, 'prefix','w')};
    
            ppparams.wcbfmap = spm_file(ppparams.qcbfmap, 'prefix','w');
        end
    end
end

%% Smoothing of the CBF map
if params.asl.do_smoothing
    if ~isfield(ppparams,'scbfmap') && (isfield(ppparams,'wcbfmap') || isfield(ppparams,'qcbfmap'))
        fprintf('Start smoothing \n')
    
        if isfield(ppparams,'wcbfmap'), Vcbf = spm_vol(ppparams.wcbfmap); else Vcbf = spm_vol(ppparams.qcbfmap); end
    
        wcbfdat = spm_read_vols(Vcbf);
    
        spm_progress_bar('Init',numel(Vcbf),'Smoothing','volumes completed');
    
        for i=1:numel(Vcbf)
            [pth,nm,~] = fileparts(Vcbf(i).fname);
    
            Q = fullfile(pth, ['s' nm  '.nii,' num2str(Vcbf(i).n)]);
            my_spmbatch_smooth(wcbfdat(:,:,:,i),Vcbf(i),Q,[params.asl.smoothfwhm params.asl.smoothfwhm params.asl.smoothfwhm],0);
    
            spm_progress_bar('Set',i);
        end
    
        spm_progress_bar('Clear');
    
        if isfield(ppparams,'wcbfmap'), ppparams.scbfmap = spm_file(ppparams.wcbfmap, 'prefix','s'); else ppparams.scbfmap = spm_file(ppparams.ecbfmap, 'prefix','s'); end
        keepfiles{numel(keepfiles)+1} = {ppparams.scbfmap};    
    
        fprintf('Done smoothing \n')
    end
end