function dirlist = my_spmbatch_dirfilelist(directory,fextension,namefilters,noprefix)

dirlist = dir(directory); %Make list of alll files
if isempty(dirlist), return; end

tmp = find(strlength({dirlist.name})>4); %Remove '.' and '..'
if ~isempty(tmp)
    dirlist = dirlist(tmp);
else
    dirlist = {};
    return 
end

tmp = find(~contains({dirlist.name},'._')); %Remove the hiden files from Mac from the list
if ~isempty(tmp)
    dirlist = dirlist(tmp);
else
    dirlist = {};
    return 
end

tmp = find(~contains({dirlist.name},'_rf_')); %Remove undeleted files from pepolar geometric correction
if ~isempty(tmp)
    dirlist = dirlist(tmp);
else
    dirlist = {};
    return 
end
tmp = find(~contains({dirlist.name},'_e_')); %Remove undeleted files from pepolar geometric correction
if ~isempty(tmp)
    dirlist = dirlist(tmp);
else
    dirlist = {};
    return 
end

if ~contains(fextension,'.'), fextension=['.' fextension]; end
tmp = find(contains({dirlist.name},fextension));%Filter list based on the file extension
if ~isempty(tmp) 
    dirlist = dirlist(tmp); 
else 
    dirlist = {};
    return 
end

if noprefix %select only those that start with sub-....
    prefixlist = split({dirlist.name},'sub-');
    prefixlist = prefixlist(:,:,1);

    tmp = find(strlength(prefixlist)==0);
    if ~isempty(tmp)
        dirlist = dirlist(tmp);
    else
        dirlist = {};
        return 
    end
end

for ifilt=1:numel(namefilters)
    if ~isempty(dirlist)
        tmp = find(contains({dirlist.name},namefilters(ifilt).name)); %Select only the files with namefilters in name
        if namefilters(ifilt).required
            if ~isempty(tmp)
                dirlist = dirlist(tmp); 
            else 
                dirlist = {};
                return 
            end
        end
    end
end