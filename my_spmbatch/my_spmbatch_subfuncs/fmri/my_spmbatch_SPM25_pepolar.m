function [ppparams,delfiles,keepfiles] = my_spmbatch_SPM25_pepolar(numdummy,ne,nt,ppparams,delfiles,keepfiles)

vdmfile = fullfile(ppparams.subfuncdir,['vdm5' ppparams.func(ne).tprefix ppparams.func(ne).funcfile]);
if exist(vdmfile,'file')>0
    ppparams.pepolar.vdm = vdmfile;
elseif nt==1 
    [Vppfunc,ppfuncdat] = my_spmbatch_readSEfMRI(ppparams.subfmapdir,ppparams.func(ne).fmapfile,numdummy+1,ppparams,1);
    
    Vppfunc.fname = fullfile(ppparams.subfmapdir,['f' ppparams.func(ne).fmapfile]);
    Vppfunc.descrip = 'my_spmbatch - first volume';
    Vppfunc.n = [1 1];
    Vppfunc = spm_create_vol(Vppfunc);
    Vppfunc = spm_write_vol(Vppfunc,ppfuncdat);
    
    auto_acpc_reorient([fullfile(ppparams.subfmapdir,['f' ppparams.func(ne).fmapfile]) ',1'],'EPI');
    MM = Vppfunc(1).mat;
    
    Vppfunc = my_reset_orientation(Vppfunc,MM);
    Vppfunc = spm_create_vol(Vppfunc);
    
    ppfunc = fullfile(ppparams.subfmapdir,['f' ppparams.func(ne).fmapfile]);
    delfiles{numel(delfiles)+1} = {ppfunc};
    
    clear ppfuncdat Vppfunc
    
    %% coregister fmap to func   
    reffunc = fullfile(ppparams.subfuncdir,[ppparams.func(ne).tprefix ppparams.func(ne).funcfile ',1']);
    
    estwrite.ref(1) = {reffunc};
    estwrite.source(1) = {ppfunc};
    estwrite.other = {ppfunc};
    estwrite.eoptions = spm_get_defaults('coreg.estimate');
    estwrite.roptions = spm_get_defaults('coreg.write');
    
    out_coreg = spm_run_coreg(estwrite);
    
    delfiles{numel(delfiles)+1} = {ppfunc};
    delfiles{numel(delfiles)+1} = {spm_file(ppfunc, 'prefix','r')};
    
    %% Using the tool SCOPE in SPM25
    fwhm = [8 4 2 1 0];
    reg = [0 10 100];
    rinterp = 2;
    jac = 1;
    prefix = 'vdm5';
    outdir = {ppparams.subfuncdir};

    vdm = spm_scope({reffunc},{out_coreg.rfiles{1}},fwhm,reg,rinterp,jac,prefix,outdir{1});
    ppparams.pepolar.vdm = vdm.dat.fname;

    delfiles{numel(delfiles)+1} = {ppparams.pepolar.vdm};
end

jsondat = fileread(ppparams.func(ne).jsonfile);
jsondat = jsondecode(jsondat);
pedir = jsondat.PhaseEncodingDirection;

if contains(pedir,'i'), ppparams.pepolar.pedim = 1;
elseif contains(pedir,'j'), ppparams.pepolar.pedim = 2;
else ppparams.pepolar.pedim = 3; end

if contains(pedir,'-')
    ppparams.pepolar.blipdir=-1;
else
    ppparams.pepolar.blipdir=1;
end