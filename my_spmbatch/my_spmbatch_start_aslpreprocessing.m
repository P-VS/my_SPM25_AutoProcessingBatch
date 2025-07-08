function my_spmbatch_start_aslpreprocessing(sublist,nsessions,datpath,params)

if ~params.asl.mruns, params.asl.runs = [1]; end

save(fullfile(datpath,'params.mat'),'params')

datlist = zeros(numel(sublist)*numel(nsessions)*numel(params.asl.runs),3);

dpos = 1;
for i = 1:numel(sublist)
    for j = 1:numel(nsessions)
        for k = 1:numel(params.asl.runs)
            datlist(dpos,1) = sublist(i);
            datlist(dpos,2) = nsessions(j);
            datlist(dpos,3) = params.asl.runs(k);
    
            dpos = dpos+1;
        end
    end
end

numpacks = ceil(numel(datlist(:,1))/params.maxprocesses);

if params.use_parallel
    for j=1:numpacks
        if (j*params.maxprocesses)<=numel(datlist(:,1))
            maxruns = params.maxprocesses;
        else
            maxruns = params.maxprocesses-((j*params.maxprocesses)-numel(datlist(:,1)));
        end

        for is = 1:maxruns
            i = (j-1)*params.maxprocesses+is;

            fprintf(['\nStart preprocessing data for subject ' num2str(datlist(i,1)) ' session ' num2str(datlist(i,2)) ' run ' num2str(datlist(i,3)) '\n'])
    
            mtlb_cmd = sprintf('"restoredefaultpath;addpath(genpath(''%s''));addpath(genpath(''%s''));my_spmbatch_run_aslpreprocessing(%d,%d,%d,''%s'',''%s'');"', ...
                                        params.spm_path,params.my_spmbatch_path,datlist(i,1),datlist(i,2),datlist(i,3),datpath,fullfile(datpath,'params.mat'));
            logfile{i} = fullfile(datpath,['asl_preprocess_logfile_' sprintf('%02d',datlist(i,1)) '_' sprintf('%02d',datlist(i,2)) '_' sprintf('%02d',datlist(i,3)) '.txt']);
    
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

                        nlogfname = fullfile(datpath,['error_asl_preprocess_logfile_' sprintf('%02d',datlist(i,1)) '_' sprintf('%02d',datlist(i,2)) '_' sprintf('%02d',datlist(i,3)) '.txt']);
                        movefile(logfile{i},nlogfname);

                        fprintf(['\nError during preprocessing data for subject ' num2str(datlist(i,1)) ' session ' num2str(datlist(i,2)) ' run ' num2str(datlist(i,3)) '\n'])
                    elseif ~isempty(test)
                        pfinnished = pfinnished+1;

                        if ~params.keeplogs
                            delete(logfile{i}); 
                        else
                            nlogfname = fullfile(datpath,['done_asl_preprocess_logfile_' sprintf('%02d',datlist(i,1)) '_' sprintf('%02d',datlist(i,2)) '_' sprintf('%02d',datlist(i,3)) '.txt']);
                            movefile(logfile{i},nlogfname);
                        end

                        fprintf(['\nDone preprocessing data for subject ' num2str(datlist(i,1)) ' session ' num2str(datlist(i,2)) ' run ' num2str(datlist(i,3)) '\n'])
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
else
    for i=1:numel(datlist(:,1))
        itstart = tic;

        my_spmbatch_run_aslpreprocessing(datlist(i,1),datlist(i,2),datlist(i,3),datpath,fullfile(datpath,'params.mat'));

        itstop = toc(itstart);

        fprintf(['subject ' num2str(datlist(i,1)) ' session ' num2str(datlist(i,2)) ' run ' num2str(datlist(i,3)) ' processed in ' datestr(duration([0,0,itstop],'InputFormat','ss'),'HH:MM:SS') '\n'])
    end
end

delete(fullfile(datpath,'params.mat'))