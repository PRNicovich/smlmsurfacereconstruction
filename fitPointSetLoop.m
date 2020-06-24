function xPiecewise = fitPointSetLoop(pts, start, stop, xSlice, axisToUse, subsampleStep, onPoint)

% Fit SMLM point cloud in slices along single plane.
% Uses graph traverse methods to approximate center
% path along a 2d point cloud. 
% Inputs:
% pts - 3D SMLM point cloud. N x 3 matrix.
% start - Start position along slicing axis. Scalar in pts units.
% stop - Stop position along slicing axis. Scalar in pts units.
% xSlice - Thickness of slice.  Scalar in pts units.
% axisToUse - Specify slice axis. 1 = slice along x axis. 2 = slice along y axis.
% subsampleStep = Smoothing factor for fitting curve.  Scalar. 
% onPoint - Number of points to skip from beginning of piecewise fitting. 
% Outputs:
% xPiecewise - cell array of N x 3 points corresponding to approximated path through point cloud in each slice.
%


xVect = start : xSlice : stop;

xPiecewise = cell(numel(xVect), 2);

switch axisToUse
    case 1
        sliceAx = 2;
        passDims = [1, 3, 4, 5];
    case 2
        sliceAx = 1;
        passDims = [2, 3, 4, 5];
        
    case 3
        sliceAx = 3;
        passDims = [1, 2, 4, 5];
end


% for xSliceNow = 84

iter = 1;

for xSliceNow = 1:numel(xVect)
    
    fprintf(1, 'Working on slice %.d of %.d.\n', xSliceNow, numel(xVect));
    
    % Pull single slice and fit
    slicePts = pts((pts(:,sliceAx) > xVect(xSliceNow)-(xSlice/2)) & (pts(:,sliceAx) < xVect(xSliceNow)+(xSlice/2)), ...
        passDims);
    
    piecewisePoints = meshfitAndSmooth2d(slicePts, subsampleStep, onPoint);

    if ~isempty(piecewisePoints)
        % piecewisePoints = meshfitAndSmooth2d(pts, 20, 5);

        % Close loop if possible
        interPointDistances = diag(squareform(pdist(piecewisePoints)), 1);

        startToEndLink = diag(squareform(pdist(piecewisePoints([1, end], :))), 1);

        % Close loop if startToEndLink is within bounds of existing
        % piecewisePoints
        % Doesn't catch every ring closure in a given object and not sure if that matters
        if ~isempty(interPointDistances)
            if (startToEndLink < max(interPointDistances(:))) && (startToEndLink > min(interPointDistances(:)))

                piecewisePoints = [piecewisePoints; piecewisePoints(1,:)];

            end
        end
    end
    
    
    xPiecewise{iter, 1} = [piecewisePoints, xVect(xSliceNow)*ones(size(piecewisePoints, 1), 1)];
    xPiecewise{iter, 2} = xVect(xSliceNow);
    
    iter = iter + 1;
    
end

