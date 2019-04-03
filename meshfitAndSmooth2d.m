function [smoothedFit, varargout] = meshfitAndSmooth2d(slicePts, subsampleStep, onPoint)

    dt = delaunayTriangulation(slicePts(:,1:2));

    connList = zeros(size(dt.ConnectivityList, 1), 2);
    
    if isempty(connList)
        smoothedFit = [];
        return
    end
    
    i = 1;
    for k = 1:size(dt.ConnectivityList, 1)

        connList(i:i+1,:) = [dt.ConnectivityList(k,1:2); ... 
                             dt.ConnectivityList(k,[1,3])];

        i = i+2;
    end

    [~, idxConn] = unique(connList, 'rows');

    startPts = dt.Points(connList(:,1),:);
    endPts = dt.Points(connList(:,2), :);

    connLength = sqrt((startPts(:,1) - endPts(:,1)).^2 + (startPts(:,2) - endPts(:,2)).^2);

    DG = sparse(connList(idxConn,1), connList(idxConn,2), connLength(idxConn), ...
        max(max(connList(idxConn, :))), max(max(connList(idxConn, :))));
    
    
%     DG = sparse(connList(:,1), connList(:,2), connLength);
    UG = tril(DG + DG');
    [ST,~] = graphminspantree(UG);

    %% Find min traverse of spanning tree

    [minA, minB, minWts] = find(ST);
    minGraph = graph(minA, minB, minWts);
    pairDist = distances(minGraph);
    [maxDistA, maxDistB] = find(pairDist == max(pairDist(:)));

    % Choose first point returned for max distance above to remove possible
    % ambiguity
    shortSpan = shortestpath(minGraph, maxDistA(1), maxDistB(1));

    minSet = dt.Points(shortSpan, :);

    %% Smooth points by built-in functions
    % Sub-sample at some reasonable spacing
    xx = smooth(1:length(minSet(:,1)), minSet(:,1), subsampleStep, 'moving');
    yy = smooth(1:length(minSet(:,2)), minSet(:,2), subsampleStep, 'moving');
% 
%     nPointsForCycle = floor((size(minSet, 1) - onPoint)/subsampleStep);

    subsampleVector = onPoint:subsampleStep:size(minSet, 1);
    piecewisePoints = [xx(subsampleVector), yy(subsampleVector)];
    
    if nargout == 1
        smoothedFit = piecewisePoints;
    elseif nargout == 2
        smoothedFit = piecewisePoints;
        varargout{1} = minSet;
    end
