function SPM_file = my_spmmbatch_tedm(SPM_file,resultmap,mask_file)

load(SPM_file);

oxM.VM = SPM.xM.VM;
oxM.xs = SPM.xM.xs;

%% === TEDM setup =====================================================

% Generate New SPM file
SPM = DefaultTEDM(SPM);

% Save parameter Info
SPM.TEDM.hist.file   = {'SPM.mat' resultmap};
SPM.TEDM.hist.prefix = '';
SPM.TEDM.hist.outfile = fullfile(resultmap,'SPM');

%% === TEDM runIADL =====================================================

%=== Sessions ===
Sess = length(SPM.nscan);

for ss = 1:Sess

    %=== Prepare Parameter ====
    c = fix(clock);
    fprintf('------------------------------------------------------------------------\n');
    fprintf('   Session %02i                                                  %2i:%2i\n',ss,c(4),c(5));
    fprintf('========================================================================\n');
    
    % Defien progress var
    Pbar = waitbar(0.0,'Setting Parameters','Name','Initialization');
    for i=1:72; fprintf('_'); end
    fprintf('\n\n');
    
    %===== MAIN PARAMETERS =====================================================
    
    %Information Progreess
    
    waitbar(0.1,Pbar,'Reading parameters','Name','Initialization');
    fprintf('   Reading parameters --------------------------- [  ]');
    
    %--- Parameters ---
    Opt.K   = SPM.TEDM.Param(ss).K;
    Opt.Del = SPM.TEDM.Param(ss).Del; 
    Opt.L_A = SPM.TEDM.Param(ss).Sp_A;
    Opt.L_F = SPM.TEDM.Param(ss).Sp_F;
    Opt.cdl = SPM.TEDM.Param(ss).cdl; 
    
    VY    = SPM.xY.VY;
    xM    = SPM.xM;
    nScan = SPM.nscan(1);
    iScan = 1; 
    
    % Check sessions
    if(ss>1)
	    % Update total number of scans
	    nScan = SPM.nscan(ss);
    
	    % Update first scan position
        iScan = 1 + sum(SPM.nscan(1:(ss-1)));    
    end
    
    DIM   = VY(iScan).dim;
    mask = spm_read_vols(spm_vol(mask_file));
    cmask = find(mask>0);
    
    fprintf('\b\b\bOk]\n');

    %===== MASCARADA ===========================================================
    waitbar(0.2,Pbar,'Get data and apply mask');
    fprintf('   Get data and apply mask ---------------------- [  ]');
    
    Vfunc = spm_vol(VY(iScan).fname);
    Y = spm_read_vols(Vfunc);
    Y = reshape(Y,[prod(DIM),numel(Vfunc)]);
    Y = Y';

    Dat = Y(:,cmask);

    % Store data
    Opt.Dat = Dat;
    
    clear('Dat','Y','Vfunc');
    
    fprintf('\b\b\bOk]\n'); 

    %=== DETRENDING =======================================================
    fprintf('   Detrending ----------------------------------- [  ]');
    waitbar(0.81,Pbar,'Detrending');
    
    AuxDat = double(Opt.Dat);
    
    % Subtract the mean
    AuxDat = tedm_SinMed(AuxDat);
    
    % Subtract trends
    AuxDat = detrend(AuxDat);
    
    % Update data
    Opt.Dat = AuxDat;
    
    clear('AuxDat');
    
    fprintf('\b\b\bOk]\n');
    waitbar(0.9,Pbar);
    
    %=== CHEKING PARAMETERS ==========================================
    fprintf('   Checking parameters -------------------------- [  ]');
    waitbar(0.95,Pbar,'Checking Parameters');
    
    %--- Defaluts --------------------------------------------
    %param.iter  = 3000;           % Number of iterations
    param.iter  = 10; %<====================================================================== Just for testing =b
    param.Ini   = 'Jdr';          % Initialization mode
    param.mgreg = 'n';            % No post-processing
    param.Preg  = 'n';            % No data reduction
    param.Verb  = 'y';            % Display verbose
    %---------------------------------------------------------
    
    param.data  = Opt.Dat;           % Data
    param.K     = Opt.K;             % Total number of components
    param.Lam   = [Opt.L_A Opt.L_F]; % Sparsity parameters
    param.Del   = Opt.Del;           % Task-related time courses
    param.cdl   = Opt.cdl;           % Similarity parameter
    
    clear('Opt');
    
    fprintf('\b\b\bOk]\n');
    for i=1:72; fprintf('_'); end
    fprintf('\n\n');
    close(Pbar);
    
    %=== CALL IADL ============================================================
    
    [D,s] = tedm_IADL(param);
    
    clear('param');
    
    %=== SAVE REGRESSORS =======================================================
    for i=1:72; fprintf('_'); end
    fprintf('\n\n');
    fprintf('   Saving results ------------------------------- [  ]');
    Pbar = waitbar(0,'Save Results...','Name','Save Results');
    
    %-Remove constant atom
    %---------------------------------------------------------------------------
    Pbar = waitbar(0,Pbar,'Remove constant atom...');
    K  = SPM.TEDM.Param(ss).K;
    iB = SPM.TEDM.Param(ss).iB;
    
    D(:,iB) = [];
    s(iB,:) = [];
    
    % Update names
    SPM.TEDM.Param(ss).names(iB) = [];
    
    K = K-1;
    
    Pbar = waitbar(1,Pbar);
    
    %-Save spatial maps
    %---------------------------------------------------------------------------
    Pbar = waitbar(0,Pbar,'Saving spatial maps...');
    
    %--- Initalise map file ---
    % Parameters, dimensions and orientation
    file_ext = '.nii';
    DIM  = SPM.xY.VY(iScan).dim;
    M    = SPM.xY.VY(iScan).mat;
    metadata = {};
    
    % Initialize spatial maps
    Vmap(1:K) = deal(struct(...
        'fname',   [],...
        'dim',     DIM,...
        'dt',      [spm_type('float32') spm_platform('bigend')],...
        'mat',     M,...
        'pinfo',   [1 0 0]',...
        'descrip', 'spm_spm:beta',...
        metadata{:}));
    
    for i = 1:K
        cmpnames        = SPM.TEDM.Param(ss).names;
        name            = [sprintf(['tedm-Ss' num2str(ss,'%02i') '_' cmpnames{i}],i) file_ext];
        Vmap(i).fname   = fullfile(resultmap,name);
        Vmap(i).descrip = sprintf('spm_spm:beta (%04d) - %s',i,name);
        
        % Save names
        SPM.TEDM.Res(ss).xS{i} = name;
    end
    
    % Write info
    Vmap = spm_data_hdr_write(Vmap);
   
    % Components
    Pbar = waitbar(0.1,Pbar);
    lns = linspace(0.1,1,K);
    
    for i = 1:K
        Cmp = NaN(size(mask));
        
        % Save Components
        Cmp(cmask) = s(i,:)';
        Vmap(i) = spm_data_write(Vmap(i),Cmp);
        
        Pbar = waitbar(lns(i),Pbar);
    end
    
    % Clear stuff
    clear('Cmp','msk','Vmap','mask','cmask');
    close(Pbar);
    
    % Clear stuff
    clear('Cmp','msk','Vmap');
    
    fprintf('\b\b\bOk]\n');
    
    %===== Update Design matrix =====
    fprintf('   Store Enhanced Design Matrix ----------------- [  ]');
    
    %--- Save Temporal Components
    SPM.TEDM.Res(ss).xD = D;
    
    %--- Selct all the regressors
    for(k = 1:K)
        SPM.TEDM.Param(ss).SetReg{k} = true;
    end
      
    fprintf('\b\b\bOk]\n');
end

%--- Create a new SPM file with the Enhacned design matrix ---

SPM = tedm_Update_fMRI_design(SPM);

SPM.xM.VM = oxM.VM;
SPM.xM.xs = oxM.xs;

save(SPM_file,"SPM")


%% === Set Defaults Settings for TEDM =====================================================

function SPM = DefaultTEDM(SPM)

% Check the number of sessions
Sess = length(SPM.nscan);

% --- General Info ---
SPM.TEDM.hist.prefix  = 'Setup_';
SPM.TEDM.hist.outfile = 'SPM';

if(Sess==1) %======== Single-Session Experiment ========

    K_A = length(SPM.xX.iC); 
    Sp_A = 85*ones(1,K_A);
    
    % Check for constrant component
    if(~isempty(SPM.xX.iB))
        
        SPM.TEDM.Param.iB = SPM.xX.iB;
        K_A  = K_A + 1;
        Sp_A = [Sp_A 0];
        
    else
        
        SPM.TEDM.Param.iB = 0;
        
    end
    
    % Update parameters
    SPM.TEDM.Touch = true;
    
    SPM.TEDM.Param.K      = K_A;
    SPM.TEDM.Param.NSrc_A = K_A;
    
    SPM.TEDM.Param.Sp_A = Sp_A;
    SPM.TEDM.Param.SpMode = 'Auto';
    
    SPM.TEDM.Param.Del  = SPM.xX.X;
    
    SPM.TEDM.Param.SimMode = 'Conservative';
    
    %--- Take names ---
    for i=1:numel(SPM.Sess.U)
        name{i} = char(SPM.Sess.U(i).name);
    end
    
    % Check constnant atom
    iB = SPM.xX.iB;
    if(~isempty(iB))
        name{iB} = 'constant';
    end
    
    SPM.TEDM.Param.names = name;

    SPM.TEDM.Param.FreeCmp = 0;
    SPM.TEDM.Param.Sp_F = [];

    SPM.TEDM.Param.Aname = SPM.TEDM.Param.names';
    
    Tcheck = repmat({true},1,SPM.TEDM.Param.NSrc_A);
    Tnames = SPM.TEDM.Param.names;
    Tdata  = num2cell(SPM.TEDM.Param.Sp_A);
    
    SPM.TEDM.Param.SpDat = [Tcheck' Tnames' Tdata'];

    SPM.TEDM.Param.cdl = tedm_AutoSimilarity(SPM,1);

else %======== Multi-session experiment =========

    % Parameters
    Dur = cumsum([0 SPM.nscan]);
    
    for ss = 1:Sess
    
        % Identify Sources
        K_A = length(SPM.Sess(ss).col); 
        Sp_A = 85*ones(1,K_A);
        
        % Check for constrant component
        if(~isempty(SPM.xX.iB(ss)))
            SPM.TEDM.Param(ss).iB = SPM.xX.iB(ss);
            K_A  = K_A + 1;
            Sp_A = [Sp_A 0];
        
        else
        
            SPM.TEDM.Param(ss).iB = 0;
        
        end
        
        % Update parameters
        SPM.TEDM.Touch(ss) = false;
        
        SPM.TEDM.Param(ss).K      = K_A;
        SPM.TEDM.Param(ss).NSrc_A = K_A;
        
        SPM.TEDM.Param(ss).Sp_A   = Sp_A;
        SPM.TEDM.Param(ss).SpMode = 'Auto';
        
        SPM.TEDM.Param(ss).SimMode = 'Conservative';
        
        % Select task-related time courses
        Cols = [SPM.Sess(ss).col, SPM.xX.iB(ss)];
        Tdur = [1+Dur(ss):Dur(ss+1)];
        
        Del = SPM.xX.X(Tdur,Cols);
        
        SPM.TEDM.Param(ss).Del = Del;
        
        % Take names
        for i=1:numel(SPM.Sess(ss).U);
            name{i} = char(SPM.Sess(ss).U(i).name);
        end
            
        % Check constant atoms
        iB = SPM.xX.iB(ss);
        if(~isempty(iB))
            name{K_A(end)} = 'constant';
        end
        
        SPM.TEDM.Param(ss).name = name;

        SPM.TEDM.Param(ss).FreeCmp = 0;
        SPM.TEDM.Param(ss).Sp_F = [];

        SPM.TEDM.Param(ss).Aname = SPM.TEDM.Param(ss).name';

        NSrc_A = handles.SPM.TEDM.Param(ss).NSrc_A;

	    Tcheck = repmat({true},1,NSrc_A);
	    Tnames = handles.SPM.TEDM.Param(ss).name;
	    Tdata  = num2cell(handles.SPM.TEDM.Param(ss).Sp_A);
    
	    SPM.TEDM.Param(ss).SpDat = [Tcheck' Tnames' Tdata'];

        SPM.TEDM.Param(ss).cdl = tedm_AutoSimilarity(SPM,ss);
    end

end