function [ppparams,delfiles,keepfiles] = my_spmbatch_preprocfunc(ppparams,params,delfiles,keepfiles)

%% Preprocessing the data per echo (reorientation,motion and geometric distortion)
for ie=ppparams.echoes
    [ppparams,delfiles,keepfiles] = my_spmbatch_preprocfunc_perecho(ppparams,params,ie,delfiles,keepfiles);

    if params.func.isaslbold && contains(ppparams.func(ie).funcfile,'_aslbold.nii') && ~contains(ppparams.func(ie).prefix,'f')
        [ppparams,delfiles,keepfiles] = my_spmbatch_split_asl_bold(params,ppparams,ie,delfiles,keepfiles);
    elseif params.func.isaslbold && contains(ppparams.func(ie).funcfile,'_aslbold.nii') && contains(ppparams.func(ie).prefix,'f')
        fname = split(ppparams.func(ie).funcfile,'_aslbold.nii');
        prefsplit = split(ppparams.func(ie).prefix,'f');
        ppparams.func(ie).perffile = ['f' prefsplit{end} fname{1} '_label.nii'];
    end
    
    if params.func.do_ArtRepair && ~contains(ppparams.func(ie).prefix,'v')
        [ppparams,delfiles,keepfiles] = my_spmbatch_artrepair(ppparams,params,ie,delfiles,keepfiles);
    end
end

ppparams.reffunc = [ppparams.func(1).prefix ppparams.func(1).funcfile ',1'];

%% Slice time correction
try
    for ie=ppparams.echoes
        if params.func.do_slicetime  && ~contains(ppparams.func(ie).prefix,'a')
            Vfunc = spm_vol(fullfile(ppparams.subfuncdir,[ppparams.func(ie).prefix ppparams.func(ie).funcfile]));
            tdim = numel(Vfunc);
            nvols = params.loadmaxvols;
    
            for ti=1:nvols:tdim
                if ti+nvols>tdim, nvols=tdim-ti+1; end
                tVfunc = Vfunc(ti:ti+nvols-1);
    
                fprintf(['Do slice time correction echo ' num2str(ie) ' vol ' num2str(ti) '-' num2str(ti+nvols-1) '\n'])
    
                funcdat = spm_read_vols(tVfunc);
                
                if ~isfield(ppparams,'SliceTimes')
                    jsondat = fileread(ppparams.func(ie).jsonfile);
                    jsondat = jsondecode(jsondat);
                
                    ppparams.tr = jsondat.RepetitionTime;
                    nsl=tVfunc(1).dim(3);
                    
                    if isfield(jsondat,'SliceTiming'), ppparams.SliceTimes = jsondat.SliceTiming; else ppparams.SliceTimes = []; end
                    if ~(numel(ppparams.SliceTimes)==nsl), ppparams.SliceTimes = []; end
    
                    if params.func.isaslbold
                        if isfield(jsondat,'LabelingDuration'), params.asl.LabelingDuration=jsondat.LabelingDuration; end
                        if isfield(jsondat,'PostLabelDelay'), params.asl.PostLabelDelay=jsondat.PostLabelDelay; end
                    end
                    
                    if isempty(ppparams.SliceTimes)
                        if isfield(jsondat,'MultibandAccelerationFactor') || isfield(jsondat,'MultibandAccellerationFactor')
                            hbf = jsondat.MultibandAccelerationFactor; 
                            nslex = ceil(nsl/hbf);
                            isl = zeros([1,nslex]);
                            isl(1:2:nslex)=[0:1:(nslex-1)/2];
                            isl(2:2:nslex)=[ceil(nslex/2):1:nslex-1];
                            isl=repmat(isl,[1,hbf]);
                            isl = isl(1:nsl);
                        else 
                            isl = [1:2:nsl 2:2:nsl];
                            isl = isl-1;
                            nslex = nsl;
                        end
                
                        if params.func.isaslbold
                            TA = ppparams.tr-params.asl.LabelingDuration-params.asl.PostLabelDelay;
                        else
                            TA = ppparams.tr;
                        end
    
                        ppparams.SliceTimes = isl*TA/nslex;
                    elseif params.func.isaslbold
                        if max(ppparams.SliceTimes)>(ppparams.tr-params.asl.LabelingDuration-params.asl.PostLabelDelay)
                            TA = ppparams.tr-params.asl.LabelingDuration-params.asl.PostLabelDelay;
        
                            ppparams.SliceTimes = ppparams.SliceTimes * TA/ppparams.tr;
                        end
                    end
    
                    if params.func.isaslbold, ppparams.SliceTimes = params.asl.LabelingDuration+params.asl.PostLabelDelay+ppparams.SliceTimes; end
                end
                
                funcdat=my_spmbatch_st(funcdat,tVfunc,ppparams.SliceTimes,ppparams.tr);

                for iv=1:nvols
                    tVfunc(iv).fname = fullfile(ppparams.subfuncdir,['a' ppparams.func(ie).prefix ppparams.func(ie).funcfile]);
                    tVfunc(iv).descrip = 'my_spmbatch';
                    tVfunc(iv).n = [ti+iv-1 1];
                
                    tVfunc(iv) = spm_write_vol(tVfunc(iv),funcdat(:,:,:,iv));
                end
            
                clear funcdat tVfunc
            end
    
            clear Vfunc
            
            ppparams.func(ie).prefix = ['a' ppparams.func(ie).prefix];
            
            delfiles{numel(delfiles)+1} = {fullfile(ppparams.subfuncdir,[ppparams.func(ie).prefix ppparams.func(ie).funcfile])};
        end
    end
catch e
    fprintf('\nPP_Error\n');
    fprintf('\nThe error was: \n%s\n',e.message)
end

%% Denoising the ME/SE fMRI data
if (params.func.isaslbold || params.func.denoise) && ~contains(ppparams.func(ppparams.echoes(1)).prefix,'d')
    
    [ppparams,delfiles,keepfiles] = my_spmbatch_fmridenoising(ppparams,params,delfiles,keepfiles);

    if params.denoise.do_DUNE
        boldfile = fullfile(ppparams.subfuncdir,[ppparams.func(params.echoes(1)).prefix ppparams.func(ppparams.echoes(1)).funcfile]);
        if ~exist(boldfile,"file"), return; end
    end
end

for ie=ppparams.echoes
    if params.func.isaslbold && contains(ppparams.func(ie).funcfile,'_aslbold.nii')
        fname = split(ppparams.func(ie).funcfile,'_aslbold.nii');
        ppparams.func(ie).funcfile = [fname{1} '_bold.nii'];
    end
end

%% Combine multiple TE timeseries for ME-fMRI
if params.func.do_echocombination && ~contains(ppparams.func(ppparams.echoes(1)).prefix,'c')
    [ppparams,delfiles] = my_spmbatch_combineMEfMRI(ppparams,params,delfiles);
end
if params.func.do_echocombination || contains(ppparams.func(ppparams.echoes(1)).prefix,'c')
    if contains(ppparams.func(1).funcfile,'_echo-')
        nfname = split(ppparams.func(1).funcfile,'_echo-');
        ppparams.func(1).funcfile = [nfname{1} '_bold.nii'];
    end
    
    ppparams.func = ppparams.func(1);
    ppparams.echoes = 1;
    ppparams.meepi = false;
end  
 
%% Normalization of the func data
for ie=ppparams.echoes
    if params.func.do_normalization && ~contains(ppparams.func(ie).prefix,'w')
        fprintf('Do normalization\n')
    
        [ppparams,delfiles,keepfiles] = my_spmbatch_normalization_func(ie,ppparams,params,delfiles,keepfiles);
    end

    %% Smooth func        
    if params.func.do_smoothing && ~contains(ppparams.func(ie).prefix,'s')
        fprintf('Do smoothing \n')

        [ppparams,delfiles,keepfiles] = my_spmbatch_dosmoothfunc(ppparams,params,ie,delfiles,keepfiles);
    end
end

for ie=ppparams.echoes
    keepfiles{numel(keepfiles)+1} = {fullfile(ppparams.subfuncdir,[ppparams.func(ie).prefix ppparams.func(ie).funcfile])}; 
end