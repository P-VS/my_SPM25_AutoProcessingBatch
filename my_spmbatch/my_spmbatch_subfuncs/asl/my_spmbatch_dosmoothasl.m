function [ppparams,delfiles,keepfiles] = my_spmbatch_dosmoothasl(ppparams,params,delfiles,keepfiles)

Vcbf = spm_vol(fullfile(ppparams.subperfdir,ppparams.perf(1).wcbffile));

tdim = numel(Vcbf);
nvols = params.loadmaxvols;
spm_progress_bar('Init',numel(Vcbf),'Smoothing','volumes completed');
for ti=1:nvols:tdim
    if ti+nvols>tdim, nvols=tdim-ti+1; end
    cbfdat = spm_read_vols(Vcbf(ti:ti+nvols-1));

    Vout = Vcbf(ti:ti+nvols-1);
    for iv=1:nvols
        scbfdat = my_spmbatch_smooth(cbfdat(:,:,:,iv),Vcbf(ti),[],[params.func.smoothfwhm params.func.smoothfwhm params.func.smoothfwhm],0);

        Vout(iv).fname = fullfile(ppparams.subperfdir,['s' ppparams.perf(1).wcbffile]);
        Vout(iv).descrip = 'my_spmbatch - smooth';
        Vout(iv).n = [ti+iv-1 1];
        Vout(iv) = spm_write_vol(Vout(iv),scbfdat);

        spm_progress_bar('Set',ti+iv-1);

        clear scbfdat
    end

    clear Vout cbfdat
end
spm_progress_bar('Clear');

keepfiles{numel(keepfiles)+1} = {fullfile(ppparams.subperfdir,['s' ppparams.perf(1).wcbffile])};    

ppparams.perf(1).scbffile = ['s' ppparams.perf(1).wcbffile];

clear Vcbf

if contains(params.asl.temp_resolution,'original')
    Vasl = spm_vol(fullfile(ppparams.subperfdir,[ppparams.perf(1).aslprefix ppparams.perf(1).aslfile]));

    tdim = numel(Vasl);
    nvols = params.loadmaxvols;
    spm_progress_bar('Init',numel(Vasl),'Smoothing','volumes completed');
    for ti=1:nvols:tdim
        if ti+nvols>tdim, nvols=tdim-ti+1; end
        asldat = spm_read_vols(Vasl(ti:ti+nvols-1));
    
        Vout = Vasl(ti:ti+nvols-1);
        for iv=1:nvols
            sasldat = my_spmbatch_smooth(asldat(:,:,:,iv),Vasl(ti),[],[params.func.smoothfwhm params.func.smoothfwhm params.func.smoothfwhm],0);
    
            Vout(iv).fname = fullfile(ppparams.subperfdir,['s' ppparams.perf(1).aslprefix ppparams.perf(1).aslfile]);
            Vout(iv).descrip = 'my_spmbatch - smooth';
            Vout(iv).n = [ti+iv-1 1];
            Vout(iv) = spm_write_vol(Vout(iv),sasldat);
    
            spm_progress_bar('Set',ti+iv-1);
    
            clear sasldat
        end
    
        clear Vout casldat
    end
    spm_progress_bar('Clear');
    
    keepfiles{numel(keepfiles)+1} = {fullfile(ppparams.subperfdir,['s' ppparams.perf(1).aslprefix ppparams.perf(1).aslfile])};    
    
    ppparams.perf(1).aslprefix = ['s' ppparams.perf(1).aslprefix];
    
    clear Vasl
end