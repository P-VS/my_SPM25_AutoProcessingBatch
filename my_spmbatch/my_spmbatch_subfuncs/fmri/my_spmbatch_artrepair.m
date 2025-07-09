function [ppparams,delfiles,keepfiles] = my_spmbatch_artrepair(ppparams,params,ne,delfiles,keepfiles)

Vfunc = spm_vol(fullfile(ppparams.subfuncdir,[ppparams.func(ne).prefix ppparams.func(ne).funcfile]));
funcdat = spm_read_vols(Vfunc);
dim = size(funcdat);

% ------------------------
% Search for bad volumes
% ------------------------ 
%  Based on Afyouni and Nichols 2018: https://www.sciencedirect.com/science/article/pii/S1053811917311229
    [DVARS,StatDVARS]=DVARSCalc(reshape(funcdat,[dim(1)*dim(2)*dim(3),dim(4)]));

    glout_idx = find(and(abs(StatDVARS.pvals)<(0.5/dim(4)),abs(StatDVARS.DeltapDvar)>5));
  
    if ~isempty(glout_idx)
        % ------------------------
        % First try to correct for rapid movement
        % ------------------------ 
        % Large intravolume motion may cause image reconstruction
        % errors, and fast motion may cause spin history effects.
        % Guess at allowable motion within a TR. For good subjects,
        % would like this value to be as low as 0.3. For clinical subjects,
        % set this threshold higher.
              mv_thresh = 1.0;  % try 0.3 for subjects with intervals with low noise
                                % try 1.0 for severely noisy subjects 
            
        mv_data = load(ppparams.rp_file);
    
        % Convert rotation movement to degrees
        mv_data(:,4:6)= mv_data(:,4:6)*180/pi; 
    
    %   % Rotation measure assumes voxel is 65 mm from origin of rotation.
        delta = zeros(dim(4),1);  % Mean square displacement in two scans
        for i = 2:dim(4)
            delta(i,1) = (mv_data(i-1,1) - mv_data(i,1))^2 +...
                    (mv_data(i-1,2) - mv_data(i,2))^2 +...
                    (mv_data(i-1,3) - mv_data(i,3))^2 +...
                    1.28*(mv_data(i-1,4) - mv_data(i,4))^2 +...
                    1.28*(mv_data(i-1,5) - mv_data(i,5))^2 +...
                    1.28*(mv_data(i-1,6) - mv_data(i,6))^2;
            delta(i,1) = sqrt(delta(i,1));
        end
        
         % Also name the scans before the big motions (v2.2 fix)
        deltaw = zeros(dim(4),1);
        for i = 1:dim(4)-1
            deltaw(i) = max(delta(i), delta(i+1));
        end
        delta(1:dim(4)-1,1) = deltaw(1:dim(4)-1,1);
      
        % Adapt the threshold  (v2.3 fix)
        delsort = sort(delta);
        if delsort(round(0.75*dim(4))) > mv_thresh
            mv_thresh = min(1.0,delsort(round(0.75*dim(4))));
        end
    
        tmvout_idx = find(delta > mv_thresh)';

        if ~isempty(tmvout_idx)
            tmp = find(ismember(tmvout_idx,glout_idx));
            if ~isempty(tmp), mvout_idx = tmvout_idx(tmp); else mvout_idx = []; end
        else 
            mvout_idx = [];
        end

        if ~isempty(mvout_idx)
            for iv=mvout_idx
                fprintf(['\nCorrect realignment vol ' num2str(iv)])
        
                Vtemp = Vfunc(iv);
                Vtemp.fname = fullfile(ppparams.subfuncdir,'temp_realignunwarp.nii');
                Vtemp.n = [1 1];
                Vtemp = spm_write_vol(Vtemp,reshape(funcdat(:,:,:,iv),[dim(1),dim(2),dim(3)]));
        
                transM = my_spmbatch_vol_set_com(Vtemp);
                transM(1:3,4) = -transM(1:3,4);
        
                Vtemp = my_reset_orientation(Vtemp,transM * Vtemp.mat);
                Vtemp = spm_create_vol(Vtemp);
        
                estwrite.ref(1) = {[Vfunc(1).fname ',1']};
                estwrite.source(1) = {[Vtemp.fname ',1']};
                estwrite.other = {[Vtemp.fname ',1']};
                estwrite.eoptions = spm_get_defaults('coreg.estimate');
                estwrite.roptions = spm_get_defaults('coreg.write');
                estwrite.eoptions.tol=estwrite.eoptions.tol*50;
                estwrite.eoptions.sep=3;
                estwrite.eoptions.cost_fun='ncc';
        
                out_coreg = spm_run_coreg(estwrite);
        
                funcdat(:,:,:,iv) = spm_read_vols(spm_vol(fullfile(ppparams.subfuncdir,'rtemp_realignunwarp.nii')));
        
                delete(fullfile(ppparams.subfuncdir,'temp_realignunwarp.nii'))
                delete(fullfile(ppparams.subfuncdir,'rtemp_realignunwarp.nii'))
                clear Vtemp out_coreg
            end

            [DVARS,StatDVARS]=DVARSCalc(reshape(funcdat,[dim(1)*dim(2)*dim(3),dim(4)]));

            glout_idx = find(and(abs(StatDVARS.pvals)<(0.5/dim(4)),abs(StatDVARS.DeltapDvar)>5));
        end
    end

    % ------------------------
    % Secondly replace bad volumes
    % ------------------------ 
    if ~isempty(glout_idx)
        jsondat = fileread(ppparams.func(ne).jsonfile);
        jsondat = jsondecode(jsondat);
        tr = jsondat.RepetitionTime;

        goodind = find(~ismember([1:dim(4)],glout_idx));
        
        if ~isempty(find(glout_idx==1))
            funcdat(:,:,:,1) = funcdat(:,:,:,goodind(1));
        end
        if ~isempty(find(glout_idx==dim(4)))
            funcdat(:,:,:,1) = funcdat(:,:,:,goodind(end));
        end

        tmp = find(and(glout_idx>1,glout_idx<dim(4)));
        if ~isempty(tmp)
            glout_idx = glout_idx(tmp);
            goodind = find(~ismember([1:dim(4)],glout_idx));

            funcdat = reshape(funcdat,[dim(1)*dim(2)*dim(3),dim(4)]);
            funcdat(:,glout_idx) = spline((goodind-1)*tr,funcdat(:,goodind),(glout_idx-1)*tr);

            funcdat = reshape(funcdat,[dim(1),dim(2),dim(3),dim(4)]);
        end
        
        V0 = Vfunc;
        for iv=1:numel(V0)
            V0(iv).fname = fullfile(ppparams.subfuncdir,['v' ppparams.func(ne).prefix ppparams.func(ne).funcfile]);
            V0(iv).descrip = 'my_spmbatch - ArtRepair';
            V0(iv).pinfo = [1,0,0];
            V0(iv).dt = [spm_type('float32'),spm_platform('bigend')];
            V0(iv).n = [iv 1];
        end

        V0 = myspm_write_vol_4d(V0,funcdat);

        delfiles{numel(delfiles)+1} = {V0.fname};
        ppparams.func(ne).prefix = ['v' ppparams.func(ne).prefix];

        clear V0
    end

clear Vfunc funcdat