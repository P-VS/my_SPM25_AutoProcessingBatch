function my_spmbatch_start_fmripreprocessing(sublist,nsessions,task,datpath,params)

params.func.isaslbold = false;
params.preprocess_asl = false;
params.asl.splitaslbold = 'none';

if ~params.func.mruns, params.func.runs = [1]; end
if params.func.meepi && numel(params.func.echoes)==1, params.func.combination='none'; end
if params.func.meepi && ~contains(params.func.combination,'none'), params.func.do_echocombination = true; else params.func.do_echocombination = false; end

if params.denoise.do_aCompCor
    params.denoise.do_noiseregression = true;
    params.denoise.do_bpfilter = true;
    params.denoise.bpfilter = [0.008 Inf];
    params.denoise.do_mot_derivatives = true;
end

if params.denoise.do_ICA_AROMA
    params.denoise.do_noiseregression = true;
    params.denoise.do_DUNE = false;
    params.denoise.do_aCompCor = true;
    params.denoise.do_bpfilter = true;
    params.denoise.bpfilter = [0.008 Inf];
    params.denoise.do_mot_derivatives = true;
end

if params.denoise.do_DUNE
    params.denoise.do_noiseregression = false;
    params.denoise.do_ICA_AROMA = false;
    params.denoise.do_aCompCor = false;
    params.denoise.do_bpfilter = true;
    params.denoise.bpfilter = [0.008 Inf];
    params.denoise.do_mot_derivatives = true;
end

save(fullfile(datpath,'params.mat'),'params')

datlist = zeros(numel(sublist)*numel(nsessions)*numel(params.func.runs),3);

dpos = 1;
for i = 1:numel(sublist)
    for j = 1:numel(nsessions)
        for k = 1:numel(params.func.runs)
            datlist(dpos,1) = sublist(i);
            datlist(dpos,2) = nsessions(j);
            datlist(dpos,3) = params.func.runs(k);
    
            dpos = dpos+1;
        end
    end
end

numpacks = ceil(numel(datlist(:,1))/params.maxprocesses);

for k = 1:numel(task)
    if params.use_parallel
        for j=1:numpacks
            if (j*params.maxprocesses)<=numel(datlist(:,1))
                maxruns = params.maxprocesses;
            else
                maxruns = params.maxprocesses-((j*params.maxprocesses)-numel(datlist(:,1)));
            end
    
            for is = 1:maxruns
                i = (j-1)*params.maxprocesses+is;
    
                fprintf(['\nStart preprocessing data for subject ' num2str(datlist(i,1)) ' session ' num2str(datlist(i,2)) ' run ' num2str(datlist(i,3)) ' task ' task{k} '\n'])
        
                mtlb_cmd = sprintf('"restoredefaultpath;addpath(genpath(''%s''));addpath(genpath(''%s''));addpath(genpath(''%s''));my_spmbatch_run_fmripreprocessing(%d,%d,%d,''%s'',''%s'',''%s'');"', ...
                                            params.GroupICAT_path,params.spm_path,params.my_spmbatch_path,datlist(i,1),datlist(i,2),datlist(i,3),task{k},datpath,fullfile(datpath,'params.mat'));
                logfile{i} = fullfile(datpath,['fmri_preprocess_logfile_' sprintf(['%0' num2str(params.sub_digits) 'd'],datlist(i,1)) '_' sprintf('%02d',datlist(i,2)) '_' sprintf('%02d',datlist(i,3)) '_' task{k} '.txt']);
        
                if exist(logfile{i},'file'), delete(logfile{i}); end
                
                if ispc
                    export_cmd = ['set PATH=' fullfile(matlabroot,'bin')];
                    [status,result] = system(export_cmd);
                    system_cmd = sprintf(['start matlab -nodesktop -nosplash -r ' mtlb_cmd ' -logfile ' logfile{i}]);
                else
                    system_cmd = sprintf([fullfile(matlabroot,'bin') '/matlab -nosplash -r ' mtlb_cmd ' -logfile ' logfile{i} ' & ']);
                end
                [status,result]=system(system_cmd);
            end
        
            %% wait for all processing to be finnished
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
    
                            nlogfname = fullfile(datpath,['error_fmri_preprocess_logfile_' sprintf(['%0' num2str(params.sub_digits) 'd'],datlist(i,1)) '_' sprintf('%02d',datlist(i,2)) '_' sprintf('%02d',datlist(i,3)) '_' task{k} '.txt']);
                            movefile(logfile{i},nlogfname);
    
                            fprintf(['\nError during preprocessing data for subject ' num2str(datlist(i,1)) ' session ' num2str(datlist(i,2)) ' run ' num2str(datlist(i,3)) ' task ' task{k} '\n'])
                        elseif ~isempty(test)
                            pfinnished = pfinnished+1;
    
                            if ~params.keeplogs
                                delete(logfile{i}); 
                            else
                                nlogfname = fullfile(datpath,['done_fmri_preprocess_logfile_' sprintf(['%0' num2str(params.sub_digits) 'd'],datlist(i,1)) '_' sprintf('%02d',datlist(i,2)) '_' sprintf('%02d',datlist(i,3)) '_' task{k} '.txt']);
                                movefile(logfile{i},nlogfname);
                            end
    
                            fprintf(['\nDone preprocessing data for subject ' num2str(datlist(i,1)) ' session ' num2str(datlist(i,2)) ' run ' num2str(datlist(i,2)) ' task ' task{k} '\n'])
                        end
                    end
                end
        
                if pfinnished==maxruns 
                    isrunning = false; 
                else
                    pause(60);
                end
            end
        end
    
        %% plot realignment parameters
        if params.preprocess_functional && params.func.do_realignment
            for i = 1:numel(datlist(:,1))
                % Print and save realignment paramers  
                save_rp_plot(datlist(i,1),datlist(i,2),datlist(i,3),task{k},datpath,params);
            end
        end
    else
        for i=1:numel(datlist(:,1))
            itstart = tic;

            my_spmbatch_run_fmripreprocessing(datlist(i,1),datlist(i,2),datlist(i,3),task{k},datpath,fullfile(datpath,'params.mat'));

            % Print and save realignment paramers  
            save_rp_plot(datlist(i,1),datlist(i,2),datlist(i,3),task{k},datpath,params);

            itstop = toc(itstart);

            fprintf(['subject ' num2str(datlist(i,1)) ' session ' num2str(datlist(i,2)) ' run ' num2str(datlist(i,3)) ' processed in ' datestr(duration([0,0,itstop],'InputFormat','ss'),'HH:MM:SS') '\n'])
        end
    end
end

delete(fullfile(datpath,'params.mat'))