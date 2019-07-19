function [xunit, yunit, zunit] = GenerateTestDataFromSTLFile(STLFileName, NPointsToGenerate, center, scaleFactor, varargin)

% Generate points randomly distributed on surface defined by STL file.
% Points are evenly distributed over surface area of object. 
%
% Inputs:
% STLFileName - path to STL file. String.
% NPointsToGenerate - number of points to scatter on surface of this STL file. Scalar.
% center - specify center of mass of generated point cloud.  1 x 3 vector.
% scaleFactor - scale STL file dimensions.  Scalar.
% 
% Output:
% [xunit, yunit, zunit] = NPointsToGenerate x 3 vector of points around surface of STL file object.
% 
% PRN 2019
% Allen Institute for Brain Science and UNSW Single Molecule Science

% STLFileName = 'D:\Dropbox\Proposals\Discovery - 3D PALM\Figures\ResolutionTestSTL.stl';

%%

if nargin == 5
    CombineCoverslipAndObjects = 1;
    coverslipFile = varargin{1};
else
    CombineCoverslipAndObjects = 0;
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Import 3D surfaces made in CAD, exported as STL file

% Enter here if coverslip and objects modeled as separate STL files
if CombineCoverslipAndObjects == 1
    
    objFV = stlread(STLFileName);
    cslipFV = stlread(coverslipFile);

    figure(2)
    patch(objFV,'FaceColor',       [0.6 .8 1.0], ...
        'EdgeColor',       'none',        ...
        'FaceLighting',    'gouraud',     ...
        'AmbientStrength', 0.15);
    hold on
    patch(cslipFV,'FaceColor',       [1 0.8 0.6], ...
        'EdgeColor',       'none',        ...
        'FaceLighting',    'gouraud',     ...
        'AmbientStrength', 0.15);
    hold off
    
elseif CombineCoverslipAndObjects == 0

    % With combined object + coverslip STL file
    objFV = stlread(STLFileName);

    figure(2)
    clf
    patch(objFV,'FaceColor',       [0.6 .8 1.0], ...
        'EdgeColor',       'none',        ...
        'FaceLighting',    'gouraud',     ...
        'AmbientStrength', 0.15);

end

% Add a camera light, and tone down the specular highlighting
camlight('headlight');
material('dull');

% Fix the axes scaling, and set a nice view angle
axis('image');
view([12 12 10]);
lightHand = findobj('parent', gca, 'type', 'light');
for k = 1:numel(lightHand)
    lightangle(lightHand(k), 110, 30)
end

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate test STORM points
% Sample all faces with likelihood equal to their fraction of the total
% area.  


% Enter here if coverslip and objects modeled as separate STL files
if CombineCoverslipAndObjects == 1
    
    % Make top of coverslip object (y value of 0) as bottom of objects and
    % across whole image area. 
    % Any faces in object that were on the coverslip cannot be selected to
    % have a molecule simulated on them
    topPoints = find(cslipFV.vertices(:,2) == 0);
    cslipTop.faces = cslipFV.faces(all(ismember(cslipFV.faces, topPoints), 2), :);
    cslipTop.vertices = cslipFV.vertices;
   
    % Find all bottom-intersecting points in object
    bottPoints = find(objFV.vertices(:,2) <= 0);
    objBott.faces = objFV.faces(all(ismember(objFV.faces, bottPoints), 2), :);
    objBott.vertices = objFV.vertices;
    
%     figure(3)
%     patch(cslipTop, 'facecolor', [1 .8 .2]);
%     hold on
%     patch(objBott, 'facecolor', [.2 .2 1], 'edgecolor', 'k');
%     hold off
%     view([0 0])
    
    % Make triangulation of non-covered coverslip surface
    exposedCslip = [cslipTop.vertices(cslipTop.faces, :); objBott.vertices(objBott.faces, :)];
    exposedCslip = unique(exposedCslip, 'rows');
    
    
    figure(3)
    plot(exposedCslip(:,1), exposedCslip(:,3), 'o', 'markersize', 2);

    triRep = delaunayTriangulation(exposedCslip(:,1), exposedCslip(:,3));
    
    % Left here, unfinished, and going with simpler plan to just have a
    % single object for objects + coverslip at this point

else
    
%     % With combined object, things below the top surface of the coverslip face
%     % need to be removed
%     lowPoints = find(objFV.vertices(:,2) < 0);
%     cslipLower.faces = objFV.faces(all(ismember(objFV.faces, lowPoints), 2), :);
%     cslipLower.vertices = objFV.vertices;
%     % objFV.vertices(cslipLower.faces, :) = [];
%     objFV.faces(all(ismember(objFV.faces, lowPoints), 2), :) = [];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sample faces to make 3D SMLM data

% http://www.mathworks.com/matlabcentral/answers/93023-is-there-a-matlab-function-that-can-compute-the-area-of-my-patch
vectA = objFV.vertices(objFV.faces(:, 2), :) - objFV.vertices(objFV.faces(:, 1), :);
vectB = objFV.vertices(objFV.faces(:, 3), :) - objFV.vertices(objFV.faces(:, 1), :);
crossVect = cross(vectA, vectB, 2);
area = 1/2 * (sqrt(sum(crossVect.^2, 2)));

% Randomly select faces to hold points
imagedFaces = randsample(1:numel(area), NPointsToGenerate, true, area);
[faceA, faceB] = histc(imagedFaces, 1:1:numel(area));
uniqueFaces = [sort(unique(imagedFaces))', faceA(faceA > 0)'];

% Within each face, choose a random position for the molecule to be
MolCoords = zeros(NPointsToGenerate, 3);
molCount = 1;

for k = 1:size(uniqueFaces, 1);
    Face = uniqueFaces(k,1);
    TriA = objFV.vertices(objFV.faces(Face, 1), :);
    TriB = objFV.vertices(objFV.faces(Face, 2), :);
    TriC = objFV.vertices(objFV.faces(Face, 3), :);
    randPoints = rand(uniqueFaces(k,2), 2);
    
    % Calculate random point within triangle patch defined by
    % objFV.vertices(objFV.faces(Face, N)) giving coordinates ABC
    % P=(1?r1???)A+(r1???(1?r2))B+(r2r1???)C
    
    samplePoints = (1 - sqrt(randPoints(:,1)))*TriA + (sqrt(randPoints(:,1)).*(1 - randPoints(:,2)))*TriB + (randPoints(:,2).*sqrt(randPoints(:,1)))*TriC;
    
    MolCoords(molCount:molCount+uniqueFaces(k,2)-1, :) = samplePoints;
    
    molCount = molCount + uniqueFaces(k,2);
    
end

% Change dimensions of MolCoords into nm, get dimension order right
MolCoords = MolCoords(:,[1 2 3])*1000;
% MolCoords(:,1:2) = MolCoords(:,1:2) + 12500;

% xunit = (MolCoords(:,1) - 0)*scaleFactor + center(1);
% yunit = (MolCoords(:,2) - 0)*scaleFactor + center(2);
% zunit = (MolCoords(:,3) - 0)*scaleFactor + center(3);

display(center)

xunit = (MolCoords(:,1) - mean(MolCoords(:,1)))*scaleFactor + center(1);
yunit = (MolCoords(:,2) - mean(MolCoords(:,2)))*scaleFactor + center(2);
zunit = (MolCoords(:,3) - mean(MolCoords(:,3)))*scaleFactor + center(3);





