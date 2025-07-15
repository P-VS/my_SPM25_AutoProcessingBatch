function my_spmbatch_start_fmriprocessing(sublist,nsessions,datpath,params)

if ~params.func.meepi, params.func.echoes = [1]; end
if ~params.func.mruns, params.func.runs = [1]; end
if ~contains(params.modality,'fmri'), params.func.echoes = [1]; end
if ~contains(params.modality,'fasl'), params.reduced_temporal_resolution = false; end
if params.reduced_temporal_resolution
    params.optimize_HRF=false; 
    params.add_derivatives = false;
    params.analysis_type = 'within_subject';
end

params.use_echoes_as_sessions = false;
if params.func.meepi && ~contains(params.fmri_prefix,'c'), params.use_echoes_as_sessions = true; end
if params.func.meepi && contains(params.fmri_prefix,'c'), params.func.echoes = [1]; end

if ~contains(params.fmri_prefix,'s'), params.do_smoothing = true; else params.do_smoothing = false; end
if contains(params.fmri_prefix,'d'), params.add_regressors = false; end

if params.func.mruns
    switch params.func.use_runs
        case 'separately'
            params.oruns = params.func.runs; 
            params.iruns = [1];
        case 'together'
            params.oruns = [1];
            params.iruns = params.func.runs;
    end
else
    params.oruns = [1];
    params.iruns = [1];
end

if params.optimize_HRF; params.add_derivatives = false; end

save(fullfile(datpath,'params.mat'),'params')

datlist = zeros(numel(sublist)*numel(nsessions),3);

dpos = 1;
for i = 1:numel(sublist)
    for j = 1:numel(nsessions)
        for k = 1:numel(params.oruns)
            datlist(dpos,1) = sublist(i);
            datlist(dpos,2) = nsessions(j);
            datlist(dpos,3) = params.oruns(k);
    
            dpos = dpos+1;
        end
    end
end

numpacks = ceil(numel(datlist(:,1))/params.maxprocesses);

for k = 1:numel(params.task)
    if params.use_parallel
        for j=1:numpacks
            if (j*params.maxprocesses)<=numel(datlist(:,1))
                maxruns = params.maxprocesses;
            else
                maxruns = params.maxprocesses-((j*params.maxprocesses)-numel(datlist(:,1)));
            end
    
            for is = 1:maxruns
                i = (j-1)*params.maxprocesses+is;
    
                t = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss');
                fprintf(['\n'  datestr(t) ' : Start processing data for subject ' num2str(datlist(i,1)) ' session ' num2str(datlist(i,2)) ' run ' num2str(datlist(i,3)) ' task ' params.task{k} '\n'])
        
                logfile{i} = fullfile(datpath,['fmri_process_logfile_' sprintf(['%0' num2str(params.sub_digits) 'd'],datlist(i,1)) '_' sprintf('%02d',datlist(i,2)) '_' sprintf('%02d',datlist(i,3)) '_' params.task{k} '.txt']);
        
                if exist(logfile{i},'file'), delete(logfile{i}); end
                
                if ispc
                    mtlb_cmd = sprintf("restoredefaultpath;addpath(genpath('%s'));addpath(genpath('%s'));my_spmbatch_run_fmriprocessing(%d,%d,%d,'%s','%s','%s');exit", ...
                                            params.spm_path,params.my_spmbatch_path,datlist(i,1),datlist(i,2),datlist(i,3),params.task{k},datpath,fullfile(datpath,'params.mat'));

                    system_cmd = sprintf(['start matlab -nodesktop -nosplash -r "%s" -logfile %s'],mtlb_cmd,logfile{i});
                else
                    mtlb_cmd = sprintf('"restoredefaultpath;addpath(genpath(''%s''));addpath(genpath(''%s''));my_spmbatch_run_fmriprocessing(%d,%d,%d,''%s'',''%s'',''%s'');exit"', ...
                                            params.spm_path,params.my_spmbatch_path,datlist(i,1),datlist(i,2),datlist(i,3),params.task{k},datpath,fullfile(datpath,'params.mat'));

                    system_cmd = sprintf([fullfile(matlabroot,'bin') '/matlab -nosplash -r ' mtlb_cmd ' -logfile ' logfile{i} ' & ']);
                end
                [status,result]=system(system_cmd);
            end
            
            %% wait for all processing to be finsished
            isrunning = true;
            pfinnished = 0;
            while isrunning
                for is = 1:maxruns
                    i = (j-1)*params.maxprocesses+is;
    
                    if exist(logfile{i},'file')
                        FID     = fopen(logfile{i},'r');
                        txt     = textscan(FID,'%s');
                        txt     = txt{1}; 
                        test=find(cellfun('isempty',strfind(txt,'PP_Completed'))==0,1,'first');
                        errortest=find(cellfun('isempty',strfind(txt,'PP_Error'))==0,1,'first');
                        fclose(FID);
    
                        if ~isempty(errortest)
                            pfinnished = pfinnished+1;
    
                            t = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss');
                            nlogfname = fullfile(datpath,['error_fmri_process_logfile_' sprintf(['%0' num2str(params.sub_digits) 'd'],datlist(i,1)) '_' sprintf('%02d',datlist(i,2)) '_' sprintf('%02d',datlist(i,3)) '_' params.task{k} '.txt']);
                            movefile(logfile{i},nlogfname);
    
                            fprintf(['\n'  datestr(t) ' : Error during processing data for subject ' num2str(datlist(i,1)) ' session ' num2str(datlist(i,2)) ' run ' num2str(datlist(i,3)) ' task ' params.task{k} '\n'])
                        elseif ~isempty(test)
                            pfinnished = pfinnished+1;
    
                            if ~params.keeplogs
                                delete(logfile{i}); 
                            else
                                nlogfname = fullfile(datpath,['done_fmri_process_logfile_' sprintf(['%0' num2str(params.sub_digits) 'd'],datlist(i,1)) '_' sprintf('%02d',datlist(i,2)) '_' sprintf('%02d',datlist(i,3)) '_' params.task{k} '.txt']);
                                movefile(logfile{i},nlogfname);
                            end
    
                            t = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss');
                            fprintf(['\n'  datestr(t) ' : Done processing data for subject ' num2str(datlist(i,1)) ' session ' num2str(datlist(i,2)) ' run ' num2str(datlist(i,3)) ' task ' params.task{k} '\n'])
                        end
                    end
                end
        
                if pfinnished>=maxruns 
                    isrunning = false; 
                else
                    pause(60);
                end
            end
        end
    else            
        for i=1:numel(datlist(:,1))
            my_spmbatch_run_fmriprocessing(datlist(i,1),datlist(i,2),datlist(i,3),params.task{k},datpath,fullfile(datpath,'params.mat'));
        end
    end
end

delete(fullfile(datpath,'params.mat'),'params')