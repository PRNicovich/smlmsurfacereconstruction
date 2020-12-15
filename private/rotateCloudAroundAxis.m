% Rotate point cloud around specified axis by specified angle.

function pts = rotateCloudAroundAxis(pts, angle, axis)

    switch lower(axis)
        case 'x'

            X = pts(:,1);
            Y = pts(:,2)*cos(angle) - pts(:,3)*sin(angle);
            Z = pts(:,2)*sin(angle) + pts(:,3)*cos(angle);
            
        case 'y'
            % Correct?
            X = pts(:,1)*cos(angle) - pts(:,3)*sin(angle);
            Y = pts(:,2);
            Z = pts(:,1)*sin(angle) + pts(:,3)*cos(angle);


        case 'z'
            
            X = pts(:,1)*cos(angle) - pts(:,2)*sin(angle);
            Y = pts(:,1)*sin(angle) + pts(:,2)*cos(angle);
            Z = pts(:,3);
            
        otherwise
            error('Axis not supported.')
    end
 
pts(:,1) = X;
pts(:,2) = Y;
pts(:,3) = Z;
