function [polygonsOut, varargout] = processWithMeshLab(inputFile, processingScript, meshlabPath, varargin)
% Process data using meshlabserver command line interface
% Inputs : inputFile - file(s) to process.  String or cell array of strings
%          processingScript - .mlx file with script to use. Output of
%               modifyMeshLabScript if modifying.
%          varargin - {1} optional array of keys to use in parsing log file.
%                           Empty cell array if no keys to use.
%                     {2} output file name.  input.STL if not provided or empty.
%              

%     meshlabPath = 'C:\program files\vcg\meshlab\meshlabserver.exe';

    [~, fb, fc] = fileparts(processingScript);

    callString = {'"$MESHLABPATH$"', ...
        '-l $PROPSFILE$', ...
        '-i $INPUTFILE$', ... 
        '-o $OUTPUTFILE$', ... 
        '-s $SCRIPTFILE$'};
    
    % Check if output file name is provided
    outFileSpecified = false;
    if nargin == 5
       outputFile = varargin{2};
       outFileString = sprintf('%s -m sa vq', outputFile);
       outFileSpecified = true;
    end
    
    % Support multiple input files
    % Input as entries in a cell array of strings
    if iscell(inputFile)
       inF = '';
       for k = 1:length(inputFile)
           if k > 1
               inF = sprintf('%s -i', inF);
           end
           inF = sprintf('%s "%s"', inF, inputFile{k});
       end
       
       if ~outFileSpecified
        outputFile = strcat(inputFile{1}(1:(end-4)), '.stl');
       end
       propsFile = strcat(inputFile{1}(1:(end-4)), '.txt');
       
    else
        inF = inputFile;
        if ~outFileSpecified
            outputFile = strcat(inputFile(1:(end-4)), '.stl');
        end
        propsFile = strcat(inputFile(1:(end-4)), '.txt');
    end

    executeString = strjoin(callString, ' ');
    executeString = strrep(executeString, '$MESHLABPATH$', meshlabPath);
    executeString = strrep(executeString, '$PROPSFILE$', propsFile);
    executeString = strrep(executeString, '$INPUTFILE$', inF);
    if outFileSpecified
        executeString = strrep(executeString, '$OUTPUTFILE$', outFileString);
    else
        executeString = strrep(executeString, '$OUTPUTFILE$', outputFile);
    end
    executeString = strrep(executeString, '$SCRIPTFILE$', strcat(fb, fc));

    display(executeString);
    
    [~,~] = dos(executeString, '-echo');
    pause(1);
    % Read in temp file

    fprintf('Reading in file %s\n', outputFile);

    if endsWith(outputFile, '.stl')
        polygonsOut = stlread_legacy(outputFile);
    elseif endsWith(outputFile, '.ply')
        polygonsOut = readSavedMeshReconstruction(outputFile);
    else
        error('Output file format not supported.');
    end

    fprintf('Parsing propsFile %s\n', propsFile);
    % Read and parse propsFile
    if nargin >= 4
        if isempty(varargin{1})
            propsOut = [];
        else
            propsOut = parsePropsFile(propsFile, varargin{1});
        end
    else
        propsOut = parsePropsFile(propsFile);
    end

    if nargout == 2
        varargout{1} = outputFile;
    elseif nargout == 3
        varargout{1} = outputFile;
        varargout{2} = propsOut;
    end
end


% Read and parse propsFile
function propsOut = parsePropsFile(propsFile, varargin)
    pID = fopen(propsFile);
    text = textscan(pID, '%s', 'delimiter', '\n');
    
    display(nargin);
    assignin('base', 'pFin', varargin);
    
    if nargin > 1
        % second input is cell array of keys to use
        if ~isempty(varargin)
            % Keys to use. Loop over them.
            for k = 1:length(varargin{1})
                key = varargin{1}{k};
                propsOut.(strrep(key, ' ', '')) = parseMeshlabLine(text, key);
            end
        else
            % No keys + don't parse anything
            fprintf(1, 'Empty keys to parse in properties file.\n');
            propsOut = [];
        end
    elseif nargin == 1
        % Use default key    
        % Look for line that starts with 'Mesh Volume'
        fprintf('Using default key in properties file.');
        key = 'Mesh Volume';
    
        propsOut.(strrep(key, ' ', '')) = parseMeshlabLine(text, key);
    end
    
    fclose(pID);
    
    function result = parseMeshlabLine(text, key)
       
        line = text{1}{~cellfun(@isempty, strfind(text{1}, key))};
        
        switch key
            case 'Mesh Volume'
                result = regexp(line, '\d*\.?\d*', 'Match');
                result = str2double(result{1});
        end
        
    end

end
