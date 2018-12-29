function [xunit, yunit, zunit] = sampleSpherePoints(nPoints, circRadius, circCenter)

    th = rand(nPoints, 1)*2*pi;
    phi = rand(nPoints, 1)*pi;
    xunit = circRadius * cos(th).* sin(phi) + circCenter(1);
    yunit = circRadius * sin(th) .* sin(phi) + circCenter(2);
    zunit = circRadius * cos(phi) + circCenter(3);
    
