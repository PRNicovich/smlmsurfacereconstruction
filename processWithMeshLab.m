function [polygonsOut, varargout] = processWithMeshLab(inputFile, processingScript)

meshlabPath = 'C:\program files\vcg\meshlab\meshlabserver.exe';

outputFile = strcat(inputFile(1:(end-4)), '.stl');

[~, fb, fc] = fileparts(processingScript);

callString = {'"$MESHLABPATH$"', ...
    '-i $INPUTFILE$', ... 
    '-o $OUTPUTFILE$', ... 
    '-s $SCRIPTFILE$'};

executeString = strjoin(callString, ' ');
executeString = strrep(executeString, '$MESHLABPATH$', meshlabPath);
executeString = strrep(executeString, '$INPUTFILE$', inputFile);
executeString = strrep(executeString, '$OUTPUTFILE$', outputFile);
executeString = strrep(executeString, '$SCRIPTFILE$', strcat(fb, fc));

[~,~] = dos(executeString, '-echo');
pause(1);
% Read in temp file

fprintf('Reading in file %s\n', outputFile);
% polygonsOut = readSavedMeshReconstruction('bunnyOutput20181228_02.ply');

% polygonsOut = pcread(outputFile);
polygonsOut = stlread(outputFile);

if nargout == 2
    varargout{1} = outputFile;
end
