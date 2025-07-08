function funcdat = my_spmbatch_noiseregression(funcdat,ppparams,params,regkind)

if ~params.preprocess_functional, params.denoise.before_normalization = false; end

mask = spm_read_vols(spm_vol(ppparams.fmask));

s = size(funcdat);
funcdat = reshape(funcdat(:,:,:,:),[prod(s(1:end-1)),s(end)]);

tmp=find(mask>0);
mfuncdat = funcdat(tmp,:);

if contains(regkind,'bold')
    if ~isfield(ppparams,'noiseregresssor')
        if params.denoise.do_ICA_AROMA
            if isfield(ppparams,'nboldica_file'), confounds = load(ppparams.nboldica_file); else confounds = []; end
        else
            if params.denoise.do_mot_derivatives && isfield(ppparams,'der_file')
                confounds = load(ppparams.der_file);
            elseif isfield(ppparams,'rp_file')
                confounds = load(ppparams.rp_file);
            else
                confounds = [];
            end
        
            if params.denoise.do_aCompCor && isfield(ppparams,'acc_file')
                acc_confounds = load(ppparams.acc_file);
    
                confounds = cat(2,confounds,acc_confounds);
            end
        end
    else
        confounds = ppparams.noiseregresssor;
    end
elseif contains(regkind,'fasl') 
    if isfield(ppparams,'naslica_file')
        confounds = load(ppparams.naslica_file);
    else
        confounds = [];
    end
end

[mfuncdat,~] = fmri_cleaning(mfuncdat(:,:),params.denoise.polort,[],confounds,[],'restoremean','on');

funcdat(tmp,:) = mfuncdat;

funcdat = reshape(funcdat(:,:),s);

clear mfuncdat tmp