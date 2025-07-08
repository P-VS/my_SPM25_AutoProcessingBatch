function [ppparams,keepfiles,delfiles] = fmri_ica_aroma(ppparams,params,keepfiles,delfiles)
%FMRI_ICA_AROMA Performs ICA-AROMA

% Initialize parameters and run ICA:

% Get t_r
jsondat = jsondecode(fileread(ppparams.func(ppparams.echoes(1)).jsonfile));
t_r = jsondat.("RepetitionTime");

ica_dir = fullfile(ppparams.subfuncdir, ['ica-dir_' ppparams.task]); % path to ICs computed at previous step

if ~exist(ica_dir,"dir")
    mkdir(ica_dir);
    
    if params.func.meepi
        for ie=ppparams.echoes
            Vtemp = spm_vol(fullfile(ppparams.subfuncdir,[ppparams.func(ie).prefix ppparams.func(ie).funcfile]));
            if ie==ppparams.echoes(1)
                voldim = Vtemp(1).dim;
                funcdat = zeros([voldim(1),voldim(2),voldim(3),numel(Vtemp)]);
            end
            funcdat = funcdat+spm_read_vols(Vtemp);
        end
    
        funcdat = funcdat ./ numel(ppparams.echoes);
    
        for iv=1:numel(Vtemp)
            Vtemp(iv).fname = fullfile(ica_dir,['c' ppparams.func(1).prefix ppparams.func(1).funcfile]);
            Vtemp(iv).descrip = 'my_spmbatch - combine';
            Vtemp(iv).pinfo = [1,0,0];
            Vtemp(iv).n = [iv 1];
        end

        Vtemp = myspm_write_vol_4d(Vtemp,funcdat);

        ica_source_file = fullfile(ica_dir,['c' ppparams.func(1).prefix ppparams.func(1).funcfile]);

        clear Vtemp funcdat

    else
        ica_source_file = fullfile(ppparams.subfuncdir,[ppparams.func(1).prefix ppparams.func(1).funcfile]);
    end
    
    %% ICA step
    
    curdir = pwd;

    fprintf('Start ICA \n') 

    do_ica(ica_source_file,ppparams.fmask, t_r, ppparams); 

    cd(curdir)

    delfiles{numel(delfiles)+1} = {fullfile(ppparams.subfuncdir, 'input_spatial_ica.m')};
end

%% AROMA (ICA classification)

fprintf('Start ICA classification\n')

[nonBOLDICdata,nonASLICdata] = ica_aroma_classification(ppparams, params, ppparams.fmask, ica_dir, t_r);

delfiles{numel(delfiles)+1} = {fullfile(ppparams.subfuncdir, 'headdata.nii')};

%% Denoising 

if ~isempty(nonBOLDICdata)
    ppparams.nboldica_file = spm_file(fullfile(ppparams.subfuncdir,[ppparams.func(ppparams.echoes(1)).prefix ppparams.func(ppparams.echoes(1)).funcfile]), 'prefix','nboldica_','ext','.txt');
    writematrix(nonBOLDICdata,ppparams.nboldica_file,'Delimiter','tab'); 

    delfiles{numel(delfiles)+1} = {ica_dir};
    delfiles{numel(delfiles)+1} = {ppparams.nboldica_file};
end

if params.func.isaslbold && contains(params.asl.splitaslbold,'meica')
    if ~isempty(nonASLICdata)
        ppparams.naslica_file = spm_file(fullfile(ppparams.subfuncdir,[ppparams.func(ppparams.echoes(1)).prefix ppparams.func(ppparams.echoes(1)).funcfile]), 'prefix','naslica_','ext','.txt');
        writematrix(nonASLICdata,ppparams.naslica_file,'Delimiter','tab'); 

        delfiles{numel(delfiles)+1} = {ppparams.naslica_file};
    end
end

clear nonBOLDICdata nonASLICdata