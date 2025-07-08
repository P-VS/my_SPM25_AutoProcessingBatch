function out = my_spmbatch_run_aslpreprocessing(sub,ses,run,datpath,paramsfile)

load(paramsfile)

try
    %% preprocess anatomical scans
    if params.preprocess_anatomical
        [delfiles,keepfiles] = my_spmbatch_preprocess_anat(sub,ses,datpath,params);
    
        % Clean up unnecessary files
        cleanup_intermediate_files(sub,ses,datpath,delfiles,keepfiles,params,'anat','preproc_anat');
    end
    
    %% preprocess asl scans
    if params.preprocess_pcasl
        [delfiles,keepfiles] = my_spmbatch_aslpreprocessed(sub,ses,run,datpath,params);
    
        % Clean up unnecessary files
        cleanup_intermediate_files(sub,ses,datpath,delfiles,keepfiles,params,'perf',params.save_folder);
    end
catch e
    fprintf('\nPP_Error\n');
    fprintf('\nThe error was: \n%s\n',e.message)
end

fprintf('\nPP_Completed\n');

out = 1;