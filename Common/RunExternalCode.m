function [out1, out2] = RunExternalCode(filename,params)
%
% This functions runs all non-Matlab programs and returns the results 
%
%   Input:  'filename' - the name of the file to run.
%           'params'   - the program parameters.  
%           
%   Output: 'out1', 'out2' - the outputs of the DOS command execution.
% 
%
% Tested: Matlab 7.6.0
% Created by: Yossi Michaeli, September 2009
% Edited by: -
%

global project_path;

[pahtstr,name,ext,ver] = fileparts(filename);
cd(pahtstr);

if (nargin == 1)
    command = [name,ext];
elseif (nargin == 2)
    command = [name,ext,params];
end

[out1, out2] = dos(command, '-echo');
if (findstr(out2,'Error') ~= 0)
    cd(project_path);
    ME = MException('RuntimeError:InternalError', out2);
    throw(ME);
end

cd(project_path);