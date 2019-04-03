function splitCells = splitByLargeDistance(xPiecewise, splitDistX, threshVal)

% Iteratively split point sets in cell of xPiecewise such that no
% consecutive points are separated by more than threshVal distance.  If an
% entry of input xPiecewise contains a distance between consecutive points
% greater than threshVal, it will be split into two entries and sever the
% long-distance link.
%
% Inputs:
% xPiecewise - cell array of N x 3 points corresponding to a center
%               approximation fit of a sliced point cloud.  Equivalent 
%               to xPiecewise or yPiecewise from MeshFitting3d.m, 
%               piecewiseVector from rosinThreshold. 
% splitDistX - cell array of vectors of distances between consecutive points in xPiecewise.
% threshVal - threshold value of distance to allow.  Distances above
%               threshVal will be split.
% %
% Outputs:
% splitCells - cell array of N x 3 points, but with no remaining consecutive 
%                 points separated by distances greater than threshVal.


extraCells = sum(splitDistX > threshVal);

xP = cell(extraCells + size(xPiecewise, 1), 1);

iter = 1;
for k = 1:size(xPiecewise, 1)
    
    splitVect = diag(squareform(pdist(xPiecewise{k, 1})), 1) > threshVal;
   
    needToSplit = any(splitVect);
    
    if needToSplit

        splitPts = find(splitVect);
        
        splitCells = mat2cell(xPiecewise{k,1}, diff([0; splitPts; numel(splitVect) + 1]));
        
        for m = 1:length(splitCells)
            
            xP{iter, 1} = splitCells{m};
            xP{iter, 2} = xPiecewise{k,2};
           
            iter = iter + 1;
            
        end
        
        
    else

        xP{iter, 1} = xPiecewise{k, 1};
        xP{iter, 2} = xPiecewise{k, 2};
        
        iter = iter + 1;
        
    end
    
end

splitCells = xP;