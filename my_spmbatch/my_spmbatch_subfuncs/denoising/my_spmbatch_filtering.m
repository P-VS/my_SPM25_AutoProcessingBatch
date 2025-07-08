function [ppparams,keepfiles,delfiles] = my_spmbatch_filtering(ppparams,params,keepfiles,delfiles)

jsondat = fileread(ppparams.func(ppparams.echoes(1)).jsonfile);
jsondat = jsondecode(jsondat);

tr = jsondat.RepetitionTime;

mask = spm_read_vols(spm_vol(ppparams.fmask));

for ie=ppparams.echoes
    if ~contains(ppparams.func(ie).prefix,'f')
        Vfunc = spm_vol(fullfile(ppparams.subfuncdir,[ppparams.func(ie).prefix ppparams.func(ie).funcfile]));
        funcdat = spm_read_vols(Vfunc);
    
        s = size(funcdat);
        funcdat = reshape(funcdat(:,:,:,:),[prod(s(1:end-1)),s(end)]);
        
        tmp=find(mask>0);
        mfuncdat = funcdat(tmp,:);

        mfuncdat = my_spm_batch_filt_funcdata(mfuncdat, params.denoise.bpfilter(1), params.denoise.bpfilter(2), tr);

        funcdat(tmp,:) = mfuncdat;
        
        funcdat = reshape(funcdat(:,:),s);
        
        clear mfuncdat
    
        for k=1:numel(Vfunc)
            Vfunc(k).fname = fullfile(ppparams.subfuncdir,['f' ppparams.func(ie).prefix ppparams.func(ie).funcfile]);
            Vfunc(k).descrip = 'my_spmbatch - bandpass filtering';
            Vfunc(k).pinfo = [1,0,0];
            Vfunc(k).n = [k 1];
        end
        
        Vfunc = myspm_write_vol_4d(Vfunc,funcdat);
    
        clear funcdat
    
        ppparams.func(ie).prefix = ['f' ppparams.func(ie).prefix];
    
        delfiles{numel(delfiles)+1} = {fullfile(ppparams.subfuncdir,[ppparams.func(ie).prefix ppparams.func(ie).funcfile])};
    end
end