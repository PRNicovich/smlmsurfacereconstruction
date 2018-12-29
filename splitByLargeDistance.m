function splitCells = splitByLargeDistance(xPiecewise, splitDistX, threshVal)

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