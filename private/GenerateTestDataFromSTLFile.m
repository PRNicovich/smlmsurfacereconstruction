function [xunit, yunit, zunit, objFV] = GenerateTestDataFromSTLFile(STLFileName, NPointsToGenerate)

% Generate points randomly distributed on surface defined by STL file.
% Points are evenly distributed over surface area of object. 
%
% Use transformSTL file to scale + translate an example STL file into
% proper size + position for processing.  
%
% Inputs:
% STLFileName - path to STL file. String.
% NPointsToGenerate - number of points to scatter on surface of this STL file. Scalar.
% 
% Output:
% [xunit, yunit, zunit] = NPointsToGenerate x 3 vector of points around surface of STL file object.
% 
% PRN 2019
% Allen Institute for Brain Science and UNSW Single Molecule Science

% STLFileName = 'D:\Dropbox\Proposals\Discovery - 3D PALM\Figures\ResolutionTestSTL.stl';

%%
objFV = stlread_legacy(STLFileName);


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

for k = 1:size(uniqueFaces, 1)
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

% Get dimension order right
MolCoords = MolCoords(:,[1 2 3]);

xunit = MolCoords(:,1);
yunit = MolCoords(:,2);
zunit = MolCoords(:,3);


