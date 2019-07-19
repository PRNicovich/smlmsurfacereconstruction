
function polyOut = readSavedMeshReconstruction(inputFile, varargin)

    % Save ply file in meshLab with these boxes checked:

    % vert : color, normal
    % Face : color

    % Additional parameters : Binary encoding NOT checked

    %% 
%     inputFile = 'C:\Users\Rusty Nicovich\Documents\MATLAB\meshFitting\bunnyOutput20181228_02.ply';

    % This gets coordinates and vertices.  Missing connectivity.
%     pcloud = pcread(inputFile);

    % Connectivity in same file, just appended after vertex list.
    fID = fopen(inputFile);
    % Load file
    allData = textscan(fID,'%s','Delimiter','\n');
    fclose(fID);

    plyFile.headerLines = [find(strcmp(allData{1}, 'ply')), find(strcmp(allData{1}, 'end_header'))];
    plyFile.vertexElements = str2double(allData{1}{~cellfun(@isempty, strfind(allData{1}, 'element vertex'))}(2+numel('element vertex'):end));
    plyFile.faceElements = str2double(allData{1}{~cellfun(@isempty, strfind(allData{1}, 'element face'))}(2+numel('element face'):end));
    plyFile.nVertexCols = numel(strfind(allData{1}{plyFile.headerLines(2) + 1}, sprintf(' ')));
    plyFile.nFaceCols = numel(strfind(allData{1}{plyFile.headerLines(2) + 1 + plyFile.vertexElements}, sprintf(' ')));

    plyFile.vertices = zeros(plyFile.vertexElements, plyFile.nVertexCols);

    for k = 1:plyFile.vertexElements
        plyFile.vertices(k,:) = cell2mat(textscan(allData{1}{plyFile.headerLines(2)+k}, ...
            strcat(repmat('%f\t', [1, (plyFile.nVertexCols-1)]), '%f\n')));
    end

    plyFile.faces = zeros(plyFile.faceElements, plyFile.nFaceCols);

    for k = 1:plyFile.faceElements
        plyFile.faces(k,:) = cell2mat(textscan(allData{1}{plyFile.headerLines(2)+plyFile.vertexElements + k}, ...
            strcat(repmat('%f\t', [1, (plyFile.nFaceCols-1)]), '%f\n')));
    end

    %%  Convert from RGB\alpha to hsv for color
    hsvColors = rgb2hsv(plyFile.faces(:,5:7)/255);

    %% Generate triangulation of mesh w/ curvature colors
    figure(3)
    clf(3)
    tri = triangulation(plyFile.faces(:,2:4)+1, plyFile.vertices(:,1), plyFile.vertices(:,2), plyFile.vertices(:,3));
    triPlot = trisurf(tri, 'edgecolor', 'none');
    set(gca,'CLim', [-.1 0.8]); % Set color limits of plot to 5th and 95th percentiles of angle measured
    set(triPlot,'FaceColor','flat',...
           'FaceVertexCData',hsvColors(:,1),...
           'CDataMapping','scaled');
    colormap('parula')

    %%
%     hold on
%     plot3(pts(:,1), pts(:,2), pts(:,3), '.', 'markersize', 1, 'markeredgecolor', [0.8, 0.8, 0.8]);
%     hold off
    
    polyOut = tri;

