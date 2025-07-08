function [ppparams,delfiles,keepfiles] = my_spmbatch_dosmoothfunc(ppparams,params,ie,delfiles,keepfiles)
 
Vfunc = spm_vol(fullfile(ppparams.subfuncdir,[ppparams.func(ie).prefix ppparams.func(ie).funcfile]));

tdim = numel(Vfunc);
nvols = params.loadmaxvols;
spm_progress_bar('Init',numel(Vfunc),'Smoothing','volumes completed');
for ti=1:nvols:tdim
    if ti+nvols>tdim, nvols=tdim-ti+1; end
    funcdat = spm_read_vols(Vfunc(ti:ti+nvols-1));

    Vout = Vfunc(ti:ti+nvols-1);
    for iv=1:nvols
        sfuncdat = my_spmbatch_smooth(funcdat(:,:,:,iv),Vfunc(ti),[],[params.func.smoothfwhm params.func.smoothfwhm params.func.smoothfwhm],0);

        Vout(iv).fname = fullfile(ppparams.subfuncdir,['s' ppparams.func(ie).prefix ppparams.func(ie).funcfile]);
        Vout(iv).descrip = 'my_spmbatch - smooth';
        Vout(iv).n = [ti+iv-1 1];
        Vout(iv) = spm_write_vol(Vout(iv),sfuncdat);

        spm_progress_bar('Set',ti+iv-1);

        clear sfuncdat
    end

    clear Vout funcdat
end
spm_progress_bar('Clear');

keepfiles{numel(keepfiles)+1} = {fullfile(ppparams.subfuncdir,['s' ppparams.func(ie).prefix ppparams.func(ie).funcfile])};    

ppparams.func(ie).prefix = ['s' ppparams.func(ie).prefix];