function [delfiles,keepfiles] = my_spmbatch_fasl(sub,ses,run,task,datpath,params)

delfiles = {};
keepfiles = {};

ppparams.sub = sub;
ppparams.ses = ses;
ppparams.run = run;
ppparams.task = task;
ppparams.use_parallel = params.use_parallel;
ppparams.save_intermediate_results = params.save_intermediate_results;
ppparams.error = false;

%% Search for the subject/session folder

ppparams.substring = ['sub-' num2str(sub,'%02d')];
if ~isfolder(fullfile(datpath,ppparams.substring)), ppparams.substring = ['sub-' num2str(sub,'%03d')]; end

ppparams.sesstring = ['ses-' num2str(ses,'%02d')];
if ~isfolder(fullfile(datpath,ppparams.substring,ppparams.sesstring)), ppparams.sesstring = ['ses-' num2str(ses,'%03d')]; end

ppparams.subpath = fullfile(datpath,ppparams.substring,ppparams.sesstring);

if ~isfolder(ppparams.subpath), ppparams.subpath = fullfile(datpath,ppparams.substring); end

if ~isfolder(ppparams.subpath)
    fprintf(['No data folder for subject ' num2str(sub) ' session ' num2str(ses) '\n'])
    fprintf('\nPP_Error\n');
    ppparams.error = true;
    return
end

ppparams.subfuncdir = fullfile(ppparams.subpath,'func');

if ~isfolder(ppparams.subfuncdir)
    fprintf(['No func data folder found for subject ' num2str(sub) ' session ' num2str(ses)])
    fprintf('\nPP_Error\n');
    ppparams.error = true;
    return
end

if params.func.pepolar
    ppparams.subfmapdir = fullfile(ppparams.subpath,'fmap');

    if ~isfolder(ppparams.subfmapdir)
        fprintf(['No fmap data folder found for subject ' num2str(sub) ' session ' num2str(ses)])
        fprintf('\nPP_Error\n');

        ppparams.subfmapdir = ppparams.subfuncdir;
    end
end

%% Step 1: Look for the needed fMRI data files

% Search for the data files
ppparams = my_spmbatch_checkfuncfiles(ppparams,params,true);

if ~isfield(ppparams,'func'), return; end

%% Step 2: Preprocessing of the ASL part of the data
if params.preprocess_asl
    ppparams.subperfdir = fullfile(ppparams.subpath,'perf');
    
    if ~isfolder(ppparams.subfuncdir)
        fprintf(['No perf data folder found for subject ' num2str(sub) ' session ' num2str(ses)])
        fprintf('\nPP_Error\n');
        ppparams.error = true;
        return
    end
    
    if ~ppparams.error
        % Search for the data files
        ppparams = my_spmbatch_checkperffiles(ppparams,params);
    
        if contains(params.asl.splitaslbold,'meica') && ~contains(ppparams.perf(1).aslprefix,'d')
            % Redo the preprocessing of the functional data
            [ppparams,delfiles,~] = my_spmbatch_preprocfunc(ppparams,params,delfiles,keepfiles);
    
            pparams.perf(1).aslprefix = ['d' pparams.perf(1).aslprefix];
        end
        
        % Do fALS preprocessing
        [delfiles,keepfiles] = my_spmbatch_preprocfasl(ppparams,params,delfiles,keepfiles);
    end
end
