function [rfuncdat,prefix,ppparams,keepfiles,delfiles] = my_spmbatch_realignunwarp(ne,nt,nvols,prefix,ppparams,params,keepfiles,delfiles)

Vref = spm_vol(fullfile(ppparams.subfuncdir,[ppparams.func(ne).tprefix ppparams.func(ne).funcfile ',1']));
Vfunc = spm_vol(fullfile(ppparams.subfuncdir,[ppparams.func(ne).tprefix ppparams.func(ne).funcfile]));

wrap = [0 0 0];
if params.func.pepolar, wrap(ppparams.pepolar.pedim) = 1; end

%% estimate the realignment parameters
if ne==params.func.echoes(1)
    funcdat = spm_read_vols(Vfunc(nt:nt+nvols-1));

    R1.dim = Vref.dim;
    R1.mat = Vref.mat;
    R1.Vol = spm_read_vols(Vref);
    R1.dt  = Vref.dt;
    R1.n   = [1 1];
     
    R1mask = my_spmbatch_mask(R1.Vol);
    R1mask_ind = find(R1mask>0.0);
    
    for iv=1:nvols    
        R2{1}(iv).dim = Vfunc(nt+iv-1).dim;
        R2{1}(iv).mat = Vfunc(nt+iv-1).mat;
        R2{1}(iv).Vol = funcdat(:,:,:,iv);
        R2{1}(iv).dt = Vfunc(nt+iv-1).dt;
        R2{1}(iv).n = Vfunc(nt+iv-1).n;

        R2mask = my_spmbatch_mask(R2{1}(iv).Vol);
        R2mask_ind = find(R2mask>0.0);
    
        mean_R1 = mean(R1.Vol(R1mask_ind),'all');
        mean_R2 = mean(R2{1}(iv).Vol(R2mask_ind),'all');

        R2{1}(iv).Vol = R2{1}(iv).Vol - mean_R2 + mean_R1;
    end

    eoptions = spm_get_defaults('realign.estimate');
    eoptions.quality = 0.95;
    eoptions.sep = 1.5;
    eoptions.rtm = 0;
    eoptions.fwhm = 1;
    eoptions.wrap = wrap;

    fprintf(['\nRealign vols ' num2str(nt) ' - ' num2str(nt+nvols-1)])
    R2 = my_spmbatch_realign(R2,R1,eoptions);

    for iv=1:nvols 
        % Save motion correction parameters
        tmpMCParam = spm_imatrix(R2{1}(iv).mat / R1.mat);
        if (nt+iv-1) == 1, ppparams.realign.offsetMCParam = tmpMCParam(1:6); end
        MP = tmpMCParam(1:6) - ppparams.realign.offsetMCParam;
    
        if ~isfield(ppparams,'rp_file')
            ppparams.rp_file = spm_file(fullfile(ppparams.subfuncdir,[ppparams.func(ne).tprefix ppparams.func(ne).funcfile]), 'prefix','rp_','ext','.txt');
        
            if params.func.meepi
                [rppath,rpname,~] = fileparts(ppparams.rp_file);
                nrpname = split(rpname,'_echo-');
                if contains(nrpname{2},'_aslbold'), ext='_aslbold'; else ext='_bold'; end
                nrp_file = fullfile(rppath,[nrpname{1} ext '.txt']);
        
                ppparams.rp_file = nrp_file;
            end
    
            save(ppparams.rp_file,'MP','-ascii');
        
            keepfiles{numel(keepfiles)+1} = {ppparams.rp_file};
        else
            save(ppparams.rp_file,'MP','-append','-ascii');
        end
    
        ppparams.realign.R(nt+iv-1).mat = R2{1}(iv).mat;
    end

    clear R1 R2 funcdat
end

%% Reslice the functional series
if ~params.func.pepolar
    roptions = spm_get_defaults('realign.write');
    roptions.wrap = wrap;
    roptions.which = 2;
    roptions.mean = 0;
    pref = 'r';

    P = Vfunc(nt:nt+nvols-1);
    for i=1:nvols
        P(i).mat = ppparams.realign.R(nt+i-1).mat;
    end

    % Reslice to reference image grid
    rfuncdat = my_spmbatch_reslice(P,Vref,roptions);
end

if params.func.pepolar
    uweoptions = spm_get_defaults('unwarp.estimate');
    uwroptions = spm_get_defaults('unwarp.write');
    uweflags.fwhm      = 0; %uweoptions.fwhm;
    uweflags.order     = uweoptions.basfcn;
    uweflags.regorder  = uweoptions.regorder;
    uweflags.lambda    = uweoptions.regwgt;
    uweflags.jm        = uwroptions.jm;
    uweflags.fot       = uweoptions.foe;
    
    if ~isempty(uweoptions.soe)
        cnt = 1;
        for i=1:size(uweoptions.soe,2)
            for j=i:size(uweoptions.soe,2)
                sotmat(cnt,1) = uweoptions.soe(i);
                sotmat(cnt,2) = uweoptions.soe(j);
                cnt = cnt+1;
            end
        end
    else
        sotmat = [];
    end
    uweflags.sot       = sotmat;
    uweflags.fwhm      = uweoptions.fwhm;
    uweflags.rem       = 0;
    uweflags.noi       = uweoptions.noi;
    uweflags.exp_round = 'First';
    
    uwrflags.interp    = 4;
    uwrflags.rem       = 0;
    uwrflags.wrap      = wrap;
    uwrflags.mask      = 1;
    uwrflags.which     = 2;
    uwrflags.mean      = 0;
    uwrflags.jm        = uwroptions.jm;
    pref    = 'ur';
    
    if uweflags.jm == 1
        uwrflags.udc = 2;
    else
        uwrflags.udc = 1;
    end
    
    sfP = spm_vol(ppparams.pepolar.vdm);
    P = Vfunc(nt:nt+nvols-1);
    for i=1:nvols
        P(i).mat = ppparams.realign.R(nt+i-1).mat;
    end

    %-Unwarp Estimate
    %--------------------------------------------------------------------------
    uweflags.sfP = sfP;
    uweflags.M = Vref.mat;

    ds = my_spmbatch_uw_estimate(P,Vref,uweflags);

    %-Unwarp Write - Sessions should be within subjects
    %--------------------------------------------------------------------------
    rfuncdat = my_spmbatch_uw_apply(cat(2,ds),Vref,uwrflags);
    rfuncdat(rfuncdat<0) = 0;
end

if nt==1, prefix = [pref prefix]; end

clear Vref Vfunc P
