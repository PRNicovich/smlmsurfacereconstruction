function [rosinOut, splitDist] = rosinThreshold(piecewiseVector, varargin)
%
% Use Rosin thresholding to determine threshold value for distance between
% consecutive points in supplied vectors. 
%
% Uses method from "Unimodal thresholding". Paul L.Rosin. (2001). 
% https://doi.org/10.1016/S0031-3203(00)00136-9
%
% Segments unimodal distribution with tail into main distribution and tail.
%
% Inputs:
% piecewiseVector - cell array, each with matrix of points in space.
%                     Code will measure distance between consecutive points
%                     and apply Rosin thresholding to that distribution of 
%                     distances. Is xPiecewise or yPiecewise from
%                     MeshFitting3d.m 
% 
% varargin - optional arguments. Only defined is first entry - boolean for
%                     providing plots in operation. 
% Outputs:
% rosinOut - threshold value of distance segmenting splitDist values
% splitDist - all distances between consecutive points in piecewiseVector.

p = inputParser;
defaultPlot = true;
addParameter(p, 'doPlot', defaultPlot);
parse(p, 'doPlot', varargin);

splitCheck = cell(size(piecewiseVector, 1), 1);
for k = 1:size(piecewiseVector, 1)
   
    splitCheck{k, 1} = diag(squareform(pdist(piecewiseVector{k, 1})), 1);
    
end

splitDist = vertcat(splitCheck{:});
[a, b] = hist(splitDist, round(numel(splitDist)/10));

% Segment with Rosin method for unimodal thresholding
[maxA, maxBin] = max(a(:));

lineSlope = (maxA - a(end))/(b(maxBin) - max(b));
intercept = maxA - lineSlope * b(maxBin);

% midPtY = (-1/lineSlope) * midPtX + b
% b - midPtY = -m*x
% b = -m*x + midPtY

% Find histogram point farthest from line defined by max occupancy bin and
% last occupied bin
v1 = [0, intercept, 0];
v2 = [b(maxBin), maxA, 0];

a1 = repmat(v1 - v2, [numel(b), 1]);
b1 = [b(:), a(:), zeros(numel(a), 1)] - repmat(v2, [numel(b), 1]);
d = sqrt(sum(cross(a1,b1,2).^2,2)) ./ sqrt(sum(a1.^2,2));
[srt, idx] = sort(d);

% Threshold is x bin of max distance from max occupancy to largest bin line
threshX = b(idx(max(srt(idx > maxBin)) == srt));

assignin('base', 'p', p);

if p.Results.doPlot{2}
    % Plot for confirmation
    figure(3)
    clf(3)
    plot(b, a, 'b');
    hold on
    plot(b(maxBin), maxA, 'rx');
    plot(max(b), a(end), 'rx');
    plot(b, b*lineSlope + intercept, 'k:')
    plot(threshX, a(b == threshX), 'ro')
    [xi,yi] = polyxpoly([threshX, threshX + max(b)], [a(b == threshX), a(b == threshX) + max(b)*(-1/lineSlope)], ...
        [b(1), b(end)], [b(1)*lineSlope + intercept, b(end)*lineSlope + intercept]);

    plot([threshX, xi], [a(b == threshX), yi], 'r:')

    hold off
    axis equal
end

rosinOut = threshX;