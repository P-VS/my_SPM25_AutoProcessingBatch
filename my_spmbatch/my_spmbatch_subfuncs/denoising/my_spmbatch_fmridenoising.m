function [ppparams,delfiles,keepfiles] = my_spmbatch_fmridenoising(ppparams,params,delfiles,keepfiles)

reffunc = [ppparams.func(ppparams.echoes(1)).prefix ppparams.func(ppparams.echoes(1)).funcfile ',1'];

%% Make masks
if ~isfield(ppparams,'fmask')
    fprintf('Make mask \n')

    Vmask = spm_vol(fullfile(ppparams.subfuncdir,reffunc));
    fmaskdat = spm_read_vols(Vmask);

    % Functional mask
    func_mask = my_spmbatch_mask(fmaskdat);
    
    Vfuncmask = Vmask(1);
    Vfuncmask.fname = fullfile(ppparams.subfuncdir,['fmask_' ppparams.func(ppparams.echoes(1)).prefix ppparams.func(ppparams.echoes(1)).funcfile]);
    Vfuncmask.descrip = 'funcmask';
    Vfuncmask = rmfield(Vfuncmask, 'pinfo'); %remove pixel info so that there is no scaling factor applied so the values
    spm_write_vol(Vfuncmask, func_mask);
    
    ppparams.fmask = Vfuncmask.fname;
    
    clear fmaskdat func_mask Vfuncmask Vmask
end

%% Do segmentation of func data
if params.denoise.do_aCompCor || params.denoise.do_ICA_AROMA || params.denoise.do_DUNE

    if ~isfield(ppparams,'fc1im') || ~isfield(ppparams,'fc2im') || ~isfield(ppparams,'fc3im')
        fprintf('Do segmentation \n')
        
        preproc.channel.vols = {fullfile(ppparams.subfuncdir,reffunc)};
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
        
        ppparams.fc1im = fullfile(ppparams.subfuncdir,['c1' ppparams.func(ppparams.echoes(1)).prefix ppparams.func(ppparams.echoes(1)).funcfile]);
        ppparams.fc2im = fullfile(ppparams.subfuncdir,['c2' ppparams.func(ppparams.echoes(1)).prefix ppparams.func(ppparams.echoes(1)).funcfile]);
        ppparams.fc3im = fullfile(ppparams.subfuncdir,['c3' ppparams.func(ppparams.echoes(1)).prefix ppparams.func(ppparams.echoes(1)).funcfile]);
        
        sname = split(reffunc,'.nii');
        delfiles{numel(delfiles)+1} = {fullfile(ppparams.subfuncdir,[sname{1} '._seg8.mat'])};
        delfiles{numel(delfiles)+1} = {fullfile(ppparams.subfuncdir,['y_' sname{1} '.nii'])};  
    end
end

%% Save files needed for DUNE
if ~params.func.denoise && params.denoise.do_DUNE
    if isfield(ppparams,'fmask'), keepfiles{numel(keepfiles)+1} = {ppparams.fmask}; end

    if isfield(ppparams,'fc1im'), keepfiles{numel(keepfiles)+1} = {ppparams.fc1im}; end
    if isfield(ppparams,'fc2im'), keepfiles{numel(keepfiles)+1} = {ppparams.fc2im}; end
    if isfield(ppparams,'fc3im'), keepfiles{numel(keepfiles)+1} = {ppparams.fc3im}; end
else
    if isfield(ppparams,'fmask'), delfiles{numel(delfiles)+1} = {ppparams.fmask}; end

    if isfield(ppparams,'fc1im'), delfiles{numel(delfiles)+1} = {ppparams.fc1im}; end
    if isfield(ppparams,'fc2im'), delfiles{numel(delfiles)+1} = {ppparams.fc2im}; end
    if isfield(ppparams,'fc3im'), delfiles{numel(delfiles)+1} = {ppparams.fc3im}; end
end

%% Extend 6 motion regresssors with derivatives and squared regressors (to 24 regressors)
if params.denoise.do_mot_derivatives && ~isfield(ppparams,'der_file')
    fprintf('Start motion derivatives \n')

    confounds = load(ppparams.rp_file);
    confounds = confounds(:,1:6);
    confounds = cat(2,confounds,cat(1,zeros(1,6),diff(confounds)));
    confounds = cat(2,confounds,power(confounds,2));

    sname = split(ppparams.rp_file,'rp_');
    ppparams.der_file = spm_file(ppparams.rp_file, 'prefix','der_','ext','.txt');

    writematrix(confounds,ppparams.der_file,'Delimiter','tab');

    keepfiles{numel(keepfiles)+1} = {ppparams.der_file};
end

%% Bandpass filering
if ~params.func.isaslbold && params.denoise.do_bpfilter
    fprintf('Start filtering \n')

    [ppparams,keepfiles,delfiles] = my_spmbatch_filtering(ppparams,params,keepfiles,delfiles);

    fprintf('Done filtering \n')
end

%% Denoising with aCompCor covariates
if params.denoise.do_aCompCor && ~isfield(ppparams,'acc_file')
    for ie=ppparams.echoes
        Vfunc = spm_vol(fullfile(ppparams.subfuncdir,[ppparams.func(ie).prefix ppparams.func(ie).funcfile]));
        if ie==ppparams.echoes(1)
            voldim = Vfunc(1).dim;
            funcdat = zeros([voldim(1),voldim(2),voldim(3),numel(Vfunc)]);
        end
        funcdat = funcdat+spm_read_vols(Vfunc);
    end

    funcdat = funcdat ./ numel(ppparams.echoes);

    fprintf('Start aCompCor \n')

    [ppparams,keepfiles] = my_spmbatch_acompcor(funcdat,ppparams,params,keepfiles);

    fprintf('Done aCompCor \n')

    clear funcdat Vfunc
end  

%% Denoise with DUNE
if params.denoise.do_DUNE
    fprintf('Start DUNE \n')

    [ppparams,keepfiles,delfiles] = my_spmbatch_dune(ppparams,params,keepfiles,delfiles);

    fprintf('Done DUNE \n')
end

%% Denoising with ICA-AROMA
if params.denoise.do_ICA_AROMA && ~isfield(ppparams,'nboldica_file')
    fprintf('Start ICA-AROMA \n')

    [ppparams,keepfiles,delfiles] = fmri_ica_aroma(ppparams,params,keepfiles,delfiles);

    fprintf('Done ICA-AROMA \n')
end

%% Noise regression
if params.denoise.do_noiseregression || params.denoise.do_ICA_AROMA
    for ie=ppparams.echoes   
        fprintf(['Do noise regression for echo ' num2str(ie) '\n'])

        Vfunc = spm_vol(fullfile(ppparams.subfuncdir,[ppparams.func(ie).prefix ppparams.func(ie).funcfile]));
        funcdat = spm_read_vols(Vfunc);

        if ~contains(ppparams.func(ie).prefix,'d')
            fbfuncdat = my_spmbatch_noiseregression(funcdat,ppparams,params,'bold');
            
            if params.func.isaslbold
                fname = split(ppparams.func(ie).funcfile,'_aslbold.nii');
                ppparams.func(ie).funcfile = [fname{1} '_bold.nii'];
            end
        
            for k=1:numel(Vfunc)
                Vfunc(k).fname = fullfile(ppparams.subfuncdir,['d' ppparams.func(ie).prefix ppparams.func(ie).funcfile]);
                Vfunc(k).descrip = 'my_spmbatch - denoise';
                Vfunc(k).pinfo = [1,0,0];
                Vfunc(k).n = [k 1];
            end
            
            Vfunc = myspm_write_vol_4d(Vfunc,fbfuncdat);
            
            ppparams.func(ie).prefix = ['d' ppparams.func(ie).prefix];

            delfiles{numel(delfiles)+1} = {fullfile(ppparams.subfuncdir,[ppparams.func(ie).prefix ppparams.func(ie).funcfile])};

            clear fbfuncdat
        end

        if params.func.isaslbold && contains(params.asl.splitaslbold,'meica')
            ppparams.subperfdir = fullfile(ppparams.subpath,'perf');
            if ~isfolder(ppparams.subperfdir), mkdir(ppparams.subperfdir); end

            fafuncdat = my_spmbatch_noiseregression(funcdat,ppparams,params,'fasl');

            fname = split(ppparams.func(ie).funcfile,'_bold.nii');
            fpref = split(ppparams.func(ie).prefix,'d');
        
            for k=1:numel(Vfunc)
                Vfunc(k).fname = fullfile(ppparams.subperfdir,['d' fpref{end} fname{1} '_asl.nii']);
                Vfunc(k).descrip = 'my_spmbatch - meica';
                Vfunc(k).pinfo = [1,0,0];
                Vfunc(k).n = [k 1];
            end
            
            Vfunc = myspm_write_vol_4d(Vfunc,fafuncdat);

            %delfiles{numel(delfiles)+1} = {fullfile(ppparams.subperfdir,['d' fpref{end} fname{1} '_asl.nii'])};

            clear fafuncdat fVfunc
        end

        clear funcdat Vfunc

        fprintf(['Done noise regresion ' num2str(ie) '\n'])
    end 
end