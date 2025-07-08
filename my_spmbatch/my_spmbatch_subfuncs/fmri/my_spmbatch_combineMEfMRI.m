function [ppparams,delfiles] = my_spmbatch_combineMEfMRI(ppparams,params,delfiles)

fprintf('Combine echoes\n')

%% Loading the fMRI time series and deleting dummy scans

nechoes = numel(ppparams.echoes);

for ie=1:nechoes
    jsondat = fileread(ppparams.func(ppparams.echoes(ie)).jsonfile);
    jsondat = jsondecode(jsondat);

    te(ie) = 1000.0*jsondat.EchoTime;

    tefunc{ie}.Vfunc = spm_vol(fullfile(ppparams.subfuncdir,[ppparams.func(ppparams.echoes(ie)).prefix ppparams.func(ppparams.echoes(ie)).funcfile]));
end

tdim = numel(tefunc{1}.Vfunc);
voldim = tefunc{1}.Vfunc(1).dim;

nfname = split(ppparams.func(1).funcfile,'_echo-');

spm_progress_bar('Init',tdim,'Combine TE images','volumes completed');

nvols = params.loadmaxvols;
for ti=1:nvols:tdim
    if ti+nvols>tdim, nvols=tdim-ti+1; end
    tefuncdat = zeros(voldim(1),voldim(2),voldim(3),nvols,nechoes);
    
    for ie=1:nechoes
        tefuncdat(:,:,:,:,ie) = spm_read_vols(tefunc{ie}.Vfunc(ti:ti+nvols-1));
    end

    switch params.func.combination
        case 'average'
            funcdat = sum(tefuncdat,5) ./ nechoes;
    
        case 'TE_weighted'
        
            sum_weights = sum(te,'all');

            funcdat = zeros(voldim(1),voldim(2),voldim(3),nvols);
            
            for ie=1:nechoes
                functidat = tefuncdat(:,:,:,:,ie);
                functidat = functidat .* te(ie);
                functidat = functidat ./ sum_weights;
        
                funcdat = funcdat+functidat; 

                clear functidat
            end
    
        case 'T2star_weighted'
            if ti==1
                tefdat = mean(tefuncdat,4);

                t2star = my_spmbacth_make_t2star_map(tefdat,te);

                VT2 = tefunc{1}.Vfunc(1);
                rmfield(VT2,'pinfo');
                VT2.fname = fullfile(ppparams.subfuncdir,[ppparams.func(1).prefix nfname{1} '_t2star.nii']);
                VT2.descrip = 'my_spmbatch - T2star map';
                VT2.dt = [spm_type('float32'),spm_platform('bigend')];
                VT2.n = [1 1];
                VT2 = spm_write_vol(VT2,t2star);

                t2star = reshape(t2star,[voldim(1)*voldim(2)*voldim(3),1]);
                t2star_ind = find(t2star>0);

                clear tefdat

                weights = zeros(voldim(1)*voldim(2)*voldim(3),nechoes);
            
                for ne=1:nechoes
                    weights(t2star_ind,ne) = repmat(-te(ne),numel(t2star_ind),1) ./ t2star(t2star_ind,1);
                    weights(t2star_ind,ne) = exp(weights(t2star_ind,ne));
                    weights(t2star_ind,ne) = repmat(te(ne),numel(t2star_ind),1) .* weights(t2star_ind,ne);
                end
            
                sum_weights = sum(weights,2);
                weights_mask = find(sum_weights>0);

                clear t2star_ind t2star
            end

            funcdat = zeros(voldim(1),voldim(2),voldim(3),nvols);

            for ne=1:nechoes
                functidat = reshape(tefuncdat(:,:,:,:,ne),[voldim(1)*voldim(2)*voldim(3),nvols]);
                functidat = functidat .* repmat(weights(:,ne),[1 nvols]);
                functidat(weights_mask,:) = functidat(weights_mask,:) ./ repmat(sum_weights(weights_mask),[1 nvols]);
        
                funcdat = funcdat+reshape(functidat,[voldim(1),voldim(2),voldim(3),nvols]); 

                clear functidat
            end 

        case 'dyn_T2star'
                funcdat = my_spmbacth_make_t2star_map(tefuncdat,te);
                
    end

    Vout = tefunc{1}.Vfunc(ti:ti+nvols-1);
    rmfield(Vout,'pinfo');
    
    for iv=1:nvols
        Vout(iv).fname = fullfile(ppparams.subfuncdir,['c' ppparams.func(1).prefix nfname{1} '_bold.nii']);
        Vout(iv).descrip = 'my_spmbatch - combine echoes';
        Vout(iv).dt = [spm_type('float32'),spm_platform('bigend')];
        Vout(iv).n = [ti+iv-1 1];
        Vout(iv) = spm_write_vol(Vout(iv),funcdat(:,:,:,iv));
    end

    clear tefuncdat funcdat Vfunc 

    spm_progress_bar('Set',ti);
end    
spm_progress_bar('Clear');

nfname = split(ppparams.func(1).funcfile,'_echo-');

delfiles{numel(delfiles)+1} = {fullfile(ppparams.subfuncdir,['c' ppparams.func(1).prefix nfname{1} '_bold.nii']);};

ppparams.func(1).funcfile = [nfname{1} '_bold.nii'];
ppparams.func(1).prefix = ['c' ppparams.func(1).prefix];
