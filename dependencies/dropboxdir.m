function [ dbd ] = dropboxdir
%DROPBOXDIR Returns the directory of the local dropbox directory.
%   DROPBOXDIR  displays the directory of the dropbox folder in the current
%   machine.
%
%   DBD = DROPBOXDIR returns the dropbox directory in the string DBD.
% 
%   Requires http://www.mathworks.com/matlabcentral/fileexchange/33381-jsonlab--a-toolbox-to-encode-decode-json-files   
% 
%   See also GETENV.
% 
% if ispc
%     userspecificpath = [getenv('homedrive'), getenv('homepath')];
% elseif isunix
%     userspecificpath = getenv('home');
% end
% 
% commonpath = [filesep, 'Dropbox (Learning Lab)', filesep];
% 
% dbd = [userspecificpath, commonpath];
% 
% end

jsonSubPath = '\Dropbox\info.json';
possiblePaths = {'APPDATA', 'LOCALAPPDATA'};

jsonPath = '';
pathsIdx = 1;
foundJson = 0;

while (~foundJson && ( pathsIdx <= length(possiblePaths)) )
   currDir = getenv(possiblePaths{pathsIdx});
   if (exist(fullfile(currDir, jsonSubPath), 'file') == 2)
       jsonPath = (fullfile(currDir, jsonSubPath));
       break
   else
       pathsIdx = pathsIdx+1;
   end
end

if ispc
    if strcmp(jsonPath, '')
       dbd = '';
       disp('Dropbox path was not found!')
    else
       jsonStruct = loadjson(jsonPath);
       dbd = jsonStruct.business.path;
    end
    dbd = fullfile(dbd,'Learning Lab Team Folder','Patlab protocols');
else
    dbd = '/Users/Tiago/Dropbox (Learning Lab)/PatonLab/';
end

end