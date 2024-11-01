function startmulticoreslave(multicoreDir)
%STARTMULTICORESLAVE  Start multi-core processing slave process.
%   STARTMULTICORESLAVE(DIRNAME) starts a slave process for function
%   STARTMULTICOREMASTER. The given directory DIRNAME is checked for data
%   files including which function to run and which parameters to use.
%
%   STARTMULTICORESLAVE (without input arguments) uses the directory
%   <TEMPDIR>/multicorefiles, where <TEMPDIR> is the directory returned by
%   function tempdir.
%
%		Markus Buehren
%		Last modified 10.04.2009
%
%   See also STARTMULTICOREMASTER.

debugMode    = 0;
showWarnings = 0;

if debugMode
    % activate all messages
    showWarnings = 1;
end

% parameters
firstWarnTime = 10;
startWarnTime = 10*60;
maxWarnTime   = 24*3600;
startWaitTime = 0.5;
maxWaitTime   = 5;

if debugMode
    firstWarnTime = 10;
    startWarnTime = 10;
    maxWarnTime   = 60;
    maxWaitTime   = 1;
end

persistent lastSessionDateStr

% get slave file directory name
if ~exist('multicoreDir', 'var') || isempty(multicoreDir)
    multicoreDir = fullfile(tempdir2, 'multicorefiles');
end
if ~exist(multicoreDir, 'dir')
    try
        mkdir(multicoreDir);
    catch
        error('Unable to create slave file directory %s.', multicoreDir);
    end
end

% initialize variables
lastEvalEndClock = clock;
lastWarnClock    = clock;
firstRun         = true;
curWarnTime      = firstWarnTime;
curWaitTime      = startWaitTime;

while 1
    parameterFileList = findfiles(multicoreDir, 'parameters_*.mat', 'nonrecursive');
    
    % get last file that is not a semaphore file
    parameterFileName = '';
    for fileNr = length(parameterFileList):-1:1
        if isempty(strfind(parameterFileList{fileNr}, 'semaphore'))
            parameterFileName = parameterFileList{fileNr};
            break % leave the for-loop
        end
    end
    
    if ~isempty(parameterFileName)
        if debugMode
            % get parameter file number for debug messages
            fileNr = str2double(regexptokens(parameterFileName,'parameters_\d+_(\d+)\.mat'));
            disp(sprintf('****** Slave is checking file nr %d *******', fileNr));
        end
        
        % load and delete last parameter file
        sem = setfilesemaphore(parameterFileName);
        loadSuccessful = true;
        if existfile(parameterFileName)
            % try to load the parameters
            lastwarn('');
            lasterror('reset');
            try
                load(parameterFileName, 'functionHandles', 'parameters'); %% file access %%
            catch
                loadSuccessful = false;
                if showWarnings
                    disp(sprintf('Warning: Unable to load parameter file %s.', parameterFileName));
                    lastMsg = lastwarn;
                    if ~isempty(lastMsg)
                        disp(sprintf('Warning message issued when trying to load:\n%s', lastMsg));
                    end
                    displayerrorstruct;
                end
            end
            
            % check if variables to load are existing
            if loadSuccessful && (~exist('functionHandles', 'var') || ~exist('parameters', 'var'))
                loadSuccessful = false;
                if showWarnings
                    disp(textwrap2(sprintf(['Warning: Either variable ''%s'' or ''%s''', ...
                        'or ''%s'' not existing after loading file %s.'], ...
                        'functionHandles', 'parameters', parameterFileName)));
                end
            end
            
            if debugMode
                if loadSuccessful
                    disp(sprintf('Successfully loaded parameter file nr %d.', fileNr));
                else
                    disp(sprintf('Problems loading parameter file nr %d.', fileNr));
                end
            end
            
            % remove parameter file
            deleteSuccessful = mbdelete(parameterFileName, showWarnings); %% file access %%
            if ~deleteSuccessful
                % If deletion is not successful it can happen that other slaves or
                % the master also use these parameters. To avoid this, ignore the
                % loaded parameters
                loadSuccessful = false;
                if debugMode
                    disp(sprintf('Problems deleting parameter file nr %d. It will be ignored', fileNr));
                end
            end
        else
            loadSuccessful = false;
            if debugMode
                disp('No parameter files found.');
            end
        end
        
        % remove semaphore and continue if loading was not successful
        if ~loadSuccessful
            removefilesemaphore(sem);
            continue
        end
        
        % Generate a temporary file which shows when the slave started working.
        % Using this file, the master can decide if the job timed out.
        % Still using the semaphore of the parameter file above.
        workingFile = strrep(parameterFileName, 'parameters', 'working');
        generateemptyfile(workingFile);
        if debugMode
            disp(sprintf('Working file nr %d generated.', fileNr));
        end
        
        % remove semaphore file
        removefilesemaphore(sem);
        
        % show progress info
        if firstRun
            disp(sprintf('First function evaluation (%s)', datestr(clock, 'mmm dd, HH:MM')));
            firstRun = false;
        elseif etime(clock, lastEvalEndClock) > 60
            disp(sprintf('First function evaluation after %s (%s)', ...
                formattime(etime(clock, lastEvalEndClock)), datestr(clock, 'mmm dd, HH:MM')));
        end
        
        %%%%%%%%%%%%%%%%%%%%%
        % evaluate function %
        %%%%%%%%%%%%%%%%%%%%%
        if debugMode
            disp(sprintf('Slave evaluates job nr %d.', fileNr));
            t0 = mbtime;
        end
        
        % Check if date string in parameter file name has changed. If yes, call
        % "clear functions" to ensure that the latest file versions are used,
        % no older versions in Matlab's memory.
        sessionDateStr = regexptokens(parameterFileName, 'parameters_(\d+)_\d+\.mat');
        if ~strcmp(sessionDateStr, lastSessionDateStr)
            clear functions
            
            if debugMode
                disp('New multicore session detected, "clear functions" called.');
            end
        end
        lastSessionDateStr = sessionDateStr;
        
        result = cell(size(parameters)); %#ok
        for k=1:numel(parameters)
            if iscell(parameters{k})
                result{k} = feval(functionHandles{k}, parameters{k}{:}); %#ok
            else
                result{k} = feval(functionHandles{k}, parameters{k}); %#ok
            end
        end
        if debugMode
            disp(sprintf('Slave finished job nr %d in %.2f seconds.', fileNr, mbtime - t0));
        end
        
        % Save result. Use file semaphore of the parameter file to reduce the
        % overhead.
        sem = setfilesemaphore(parameterFileName);
        resultFileName = strrep(parameterFileName, 'parameters', 'result');
        try
            save(resultFileName, 'result'); %% file access %%
            if debugMode
                disp(sprintf('Result file nr %d generated.', fileNr));
            end
        catch
            if showWarnings
                disp(sprintf('Warning: Unable to save file %s.', resultFileName));
                displayerrorstruct;
            end
        end
        
        % remove working file
        mbdelete(workingFile, showWarnings); %% file access %%
        if debugMode
            disp(sprintf('Working file nr %d deleted.', fileNr));
        end
        
        % remove parameter file (might have been re-generated again by master)
        mbdelete(parameterFileName, showWarnings); %% file access %%
        if debugMode
            disp(sprintf('Parameter file nr %d deleted.', fileNr));
        end
        
        % remove semaphore
        removefilesemaphore(sem);
        
        % save time
        lastEvalEndClock = clock;
        curWarnTime = startWarnTime;
        curWaitTime = startWaitTime;
        
        % remove variables before next run
        clear result functionHandle parameters
        
    else
        % display message if idle for long time
        timeSinceLastEvaluation = etime(clock, lastEvalEndClock);
        if min(timeSinceLastEvaluation, etime(clock, lastWarnClock)) > curWarnTime
            if timeSinceLastEvaluation >= 10*60
                % round to minutes
                timeSinceLastEvaluation = 60 * round(timeSinceLastEvaluation / 60);
            end
            disp(sprintf('Warning: No slave files found during last %s (%s).', ...
                formattime(timeSinceLastEvaluation), datestr(clock, 'mmm dd, HH:MM')));
            lastWarnClock = clock;
            if firstRun
                curWarnTime = startWarnTime;
            else
                curWarnTime = min(curWarnTime * 2, maxWarnTime);
            end
            curWaitTime = min(curWaitTime + 0.5, maxWaitTime);
        end
        
        % wait before next check
        pause(curWaitTime);
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function timeString = formattime(time, mode)
%FORMATTIME  Return formatted time string.
%   STR = FORMATTIME(TIME) returns a formatted time string for the given
%   time difference TIME in seconds, i.e. '1 hour and 5 minutes' for TIME =
%   3900.
%
%   FORMATTIME(TIME, MODE) uses the specified display mode ('long' or
%   'short'). Default is long display.
%
%   FORMATTIME (without input arguments) shows examples.

if nargin == 0
    disp(sprintf('\nExamples for strings returned by function %s.m:', mfilename));
    time = [0 1e-4 0.1 1 1.1 2 60 61 62 120 121 122 3600 3660 3720 7200 7260 7320 ...
        3600*24 3600*25 3600*26 3600*48 3600*49 3600*50];
    for k=1:length(time)
        disp(sprintf('time = %6g, timeString = ''%s''', time(k), formattime(time(k))));
    end
    if nargout > 0
        timeString = '';
    end
    return
end

if ~exist('mode', 'var')
    mode = 'long';
end

if time < 0
    disp('Warning: Time must be greater or equal zero.');
    timeString = '';
elseif time >= 3600*24
    days = floor(time / (3600*24));
    if days > 1
        dayString = 'days';
    else
        dayString = 'day';
    end
    hours = floor(mod(time, 3600*24) / 3600);
    if hours == 0
        timeString = sprintf('%d %s', days, dayString);
    else
        if hours > 1
            hourString = 'hours';
        else
            hourString = 'hour';
        end
        timeString = sprintf('%d %s and %d %s', days, dayString, hours, hourString);
    end
    
elseif time >= 3600
    hours = floor(mod(time, 3600*24) / 3600);
    if hours > 1
        hourString = 'hours';
    else
        hourString = 'hour';
    end
    minutes = floor(mod(time, 3600) / 60);
    if minutes == 0
        timeString = sprintf('%d %s', hours, hourString);
    else
        if minutes > 1
            minuteString = 'minutes';
        else
            minuteString = 'minute';
        end
        timeString = sprintf('%d %s and %d %s', hours, hourString, minutes, minuteString);
    end
    
elseif time >= 60
    minutes = floor(time / 60);
    if minutes > 1
        minuteString = 'minutes';
    else
        minuteString = 'minute';
    end
    seconds = floor(mod(time, 60));
    if seconds == 0
        timeString = sprintf('%d %s', minutes, minuteString);
    else
        if seconds > 1
            secondString = 'seconds';
        else
            secondString = 'second';
        end
        timeString = sprintf('%d %s and %d %s', minutes, minuteString, seconds, secondString);
    end
    
else
    if time > 10
        seconds = floor(time);
    else
        seconds = floor(time * 100) / 100;
    end
    if seconds > 0
        if seconds ~= 1
            timeString = sprintf('%.4g seconds', seconds);
        else
            timeString = '1 second';
        end
    else
        timeString = sprintf('%.4g seconds', time);
    end
end

switch mode
    case 'long'
        % do nothing
    case 'short'
        timeString = strrep(timeString, ' and ', ' ');
        timeString = strrep(timeString, ' days', 'd');
        timeString = strrep(timeString, ' day', 'd');
        timeString = strrep(timeString, ' hours', 'h');
        timeString = strrep(timeString, ' hour', 'h');
        timeString = strrep(timeString, ' minutes', 'm');
        timeString = strrep(timeString, ' minute', 'm');
        timeString = strrep(timeString, ' seconds', 's');
        timeString = strrep(timeString, ' second', 's');
    otherwise
        error('Mode ''%s'' unknown in function %s.', mode, mfilename);
end


