function datFile = fix_dat_path(DataPath)

    configFile = fullfile(DataPath, 'config.mat');

    if ~exist(configFile, 'file')
        error('config.mat not found in: %s', DataPath);
    end

    load(configFile, 'ops')

    parentPath = fileparts(DataPath);

    d = dir(fullfile(parentPath, '*concat_CLEAN.dat'));

    if isempty(d)
        d = dir(fullfile(parentPath, '*concat*.dat'));
    end

    if isempty(d)
        error('No concat .dat file found in: %s', parentPath);
    end

    [~, idx] = max([d.bytes]);
    datFile = fullfile(d(idx).folder, d(idx).name);

    ops.fbinary  = datFile;
    ops.fclean   = datFile;
    ops.Nchan    = ops.NchanTOT;
    ops.datatype = 'int16';

    save(configFile, 'ops')

    fprintf('Using dat file: %s\n', datFile);
end
% function datFile = fix_dat_path(DataPath)
% 
%     configFile = fullfile(DataPath, 'config.mat');
% 
%     if ~exist(configFile, 'file')
%         error('config.mat not found in: %s', DataPath);
%     end
% 
%     load(configFile, 'ops')
% 
%     parentPath = fileparts(DataPath);
% 
%     d = dir(fullfile(parentPath, '*concat_CLEAN.dat'));
% 
%     if isempty(d)
%         d = dir(fullfile(parentPath, '*concat*.dat'));
%     end
% 
%     if isempty(d)
%         error('No concat .dat file found in: %s', parentPath);
%     end
% 
%     [~, idx] = max([d.bytes]);
%     datFile = fullfile(d(idx).folder, d(idx).name);
% 
%     ops.fbinary  = datFile;
%     ops.fclean   = datFile;
%     ops.Nchan    = ops.NchanTOT;
%     ops.datatype = 'int16';
% 
%     save(configFile, 'ops')
% 
%     EPA requires a .dat physically in DataPath, so create link/copy there
%     localDat = fullfile(DataPath, d(idx).name);
% 
%     if ~exist(localDat, 'file')
%         fprintf('Creating local .dat link in kilosort folder...\n')
% 
%         Try Windows symbolic link first
%         cmd = sprintf('mklink "%s" "%s"', localDat, datFile);
%         [status, msg] = system(cmd);
% 
%         If symlink fails, fall back to copying
%         if status ~= 0
%             warning('Symlink failed, copying .dat instead. This may take a while.');
%             copyfile(datFile, localDat);
%         end
%     end
% 
%     fprintf('Using dat file: %s\n', datFile);
% end