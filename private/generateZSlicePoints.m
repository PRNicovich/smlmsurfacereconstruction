% Given set of points along X or Y axis, find all intersections of segments
% with Z axis. 

function zPout = generateZSlicePoints(xPiecewise, zPlane, orientation)

switch orientation
    case 'x'
        orMat = [3 1 2];
    case 'y'
        orMat = [1 3 2];
    otherwise
        error('Orientation not supported')
end

    zP = [];
    for x = 1:length(xPiecewise)

       if ~isempty(xPiecewise{x, 1})

           transPts = find(diff(xPiecewise{x, 1}(:,2) > zPlane) ~= 0);

           if ~isempty(transPts)

               for m = 1:numel(transPts)

                   oneSide = xPiecewise{x, 1}(transPts(m), orMat);
                   otherSide = xPiecewise{x, 1}(transPts(m)+1, orMat);

                   rat = 1 - ((oneSide(3) - zPlane)/(oneSide(3) - otherSide(3)));

                   zP = [zP; rat*(oneSide - otherSide) + otherSide];

               end

           end

       end

    end
    
zPout = zP;