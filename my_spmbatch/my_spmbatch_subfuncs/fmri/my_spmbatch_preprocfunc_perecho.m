function [ppparams,delfiles,keepfiles] = my_spmbatch_preprocfunc_perecho(ppparams,params,enum,delfiles,keepfiles)

ppparams.func(enum).tprefix = '';
prefix = '';

if ~contains(ppparams.func(enum).prefix,'e'), do_reorientation = true; else 
    do_reorientation = false; 
    ppparams.func(enum).tprefix = ['e' ppparams.func(enum).tprefix];
end
if (params.func.do_realignment && ~contains(ppparams.func(enum).prefix,'r')), do_realignment = true; else 
    do_realignment = false; 
    if params.func.do_realignment, ppparams.func(enum).tprefix = ['r' ppparams.func(enum).tprefix]; end
end
if params.func.pepolar && ~contains(ppparams.func(enum).prefix,'u'), do_pepolar = true; else 
    do_pepolar = false; 
    if params.func.pepolar, ppparams.func(enum).tprefix = ['u' ppparams.func(enum).tprefix]; end
end

jsondat = fileread(ppparams.func(enum).jsonfile);
jsondat = jsondecode(jsondat);
tr = jsondat.RepetitionTime;

numdummy = floor(params.func.dummytime/tr);

Vfunc = spm_vol(fullfile(ppparams.subfuncdir,ppparams.func(enum).funcfile));
Vfunc = Vfunc(numdummy+1:end);

tdim = numel(Vfunc);
nvols = params.loadmaxvols;
for ti=1:nvols:tdim
    if ti+nvols>tdim, nvols=tdim-ti+1; end
    fprintf(['Echo: ' num2str(enum) ' vols: ' num2str(ti) '-' num2str(ti+nvols-1) ': ' '\n'])

    %% Load func data and remove dummy scans
    if do_reorientation
        [tVfunc,funcdat] = my_spmbatch_readSEfMRI(ppparams.subfuncdir,ppparams.func(enum).funcfile,numdummy+ti,ppparams,nvols);

        dim = size(funcdat);
        if numel(dim)<4, reshape(funcdat,[dim(1),dim(2),dim(3),1]); end
         
        for iv=1:nvols
            tVfunc(iv).fname = fullfile(ppparams.subfuncdir,['e' ppparams.func(enum).funcfile]);
            tVfunc(iv).descrip = 'my_spmbatch - reorient';
            tVfunc(iv).n = [ti+iv-1 1];
            tVfunc(iv) = spm_write_vol(tVfunc(iv),funcdat(:,:,:,iv));
        end
    
        if ti==1
            delfiles{numel(delfiles)+1} = {fullfile(ppparams.subfuncdir,['e' ppparams.func(enum).funcfile])}; 
            ppparams.func(enum).tprefix = ['e' ppparams.func(enum).tprefix];
        end
            
        % set orientation according to the AC-PC line   
        if enum==params.func.echoes(1) && ti==1
            auto_acpc_reorient([fullfile(ppparams.subfuncdir,[ppparams.func(enum).tprefix ppparams.func(enum).funcfile]) ',1'],'EPI');
            t2Vfunc = spm_vol(fullfile(ppparams.subfuncdir,[ppparams.func(enum).tprefix ppparams.func(enum).funcfile ',1']));
            ppparams.MM = t2Vfunc.mat;
        elseif ~isfield(ppparams,'MM')
            t2Vfunc = spm_vol(fullfile(ppparams.subfuncdir,[ppparams.func(enum).tprefix ppparams.func(enum).funcfile ',1']));
            ppparams.MM = t2Vfunc.mat;
        end
        
        for iv=1:nvols
            tVfunc(iv) = my_reset_orientation(tVfunc(iv),ppparams.MM);
            tVfunc(iv) = spm_create_vol(tVfunc(iv));
        end

        clear funcdat t2Vfunc
    end

    %% Topup geometric correction
    if do_pepolar
        [ppparams,delfiles,keepfiles] = my_spmbatch_SPM25_pepolar(numdummy,enum,ti,ppparams,delfiles,keepfiles);
    end

    %% Realignment
    if do_realignment
        [funcdat,prefix,ppparams,keepfiles,delfiles] = my_spmbatch_realignunwarp(enum,ti,nvols,prefix,ppparams,params,keepfiles,delfiles);
    end
    
    if ti==1, spm_progress_bar('Init',tdim,'Realignment','volumes completed'); end
    
    %% Save result
    if exist('funcdat','var')
        tVfunc = spm_vol(fullfile(ppparams.subfuncdir,[ppparams.func(enum).tprefix ppparams.func(enum).funcfile]));
        mat = tVfunc(1).mat;
        
        % in ASL-BOLD the last dynamic is the M0 image 
        if (ti+nvols-1)==tdim && params.func.isaslbold && contains(params.asl.isM0scan,'last')
            if ~exist(fullfile(ppparams.subpath,'perf'),'dir'), mkdir(fullfile(ppparams.subpath,'perf')); end
    
            M0 = tVfunc(end);
            m0dat = funcdat(:,:,:,end);
    
            fparts = split(ppparams.func(enum).funcfile,'_bold.nii');
            fparts = split(fparts{1},'_aslbold.nii');
            
            M0         = rmfield(M0, {'private'});
            M0.fname   = fullfile(ppparams.subpath,'perf',[prefix ppparams.func(enum).tprefix fparts{1} '_m0scan.nii']);
            M0.descrip = 'my_spmbatch - m0scan';
            M0.n       = [1 1];
            M0.mat     = mat;
            M0         = spm_create_vol(M0);
            M0         = spm_write_vol(M0,m0dat);

            M0 = my_reset_orientation(M0,mat);

            clear M0 m0dat

            tVfunc = tVfunc(1:end-1);
            funcdat = funcdat(:,:,:,1:end-1);
            nvols = nvols-1;

            if do_realignment && enum==params.func.echoes(1)
                confounds = load(ppparams.rp_file);
                confounds = confounds(1:end-1,:);

                delete(ppparams.rp_file)

                save(ppparams.rp_file,'confounds','-ascii');
            end
        end
      
        tVfunc = tVfunc(ti:ti+nvols-1);

        dim = size(funcdat);
        if numel(dim)<4, reshape(funcdat,[dim(1),dim(2),dim(3),1]); end

        for iv=1:nvols
            tVfunc(iv).fname   = fullfile(ppparams.subfuncdir,[prefix ppparams.func(enum).tprefix ppparams.func(enum).funcfile]);
            tVfunc(iv).descrip = 'my_spmbatch';
            tVfunc(iv).mat     = mat;
            tVfunc(iv).n       = [ti+iv-1 1];
        
            tVfunc(iv) = spm_write_vol(tVfunc(iv),funcdat(:,:,:,iv));
        end

        tVfunc(iv) = my_reset_orientation(tVfunc(iv),mat);
    
        clear funcdat tVfunc
    end

    spm_progress_bar('Set',ti+nvols-1);
end
spm_progress_bar('Clear');

if do_reorientation || do_realignment || do_pepolar
    delfiles{numel(delfiles)+1} = {fullfile(ppparams.subfuncdir,[prefix ppparams.func(enum).tprefix ppparams.func(enum).funcfile])};
    ppparams.func(enum).prefix = [prefix ppparams.func(enum).tprefix];
end