
fileName = 'smlm3dreconBunny.gif';

 gifDelay = 0.05;
% Original
h = figure(2);
clf
set(h, 'position', [705.6667  291.0000  600.0000  400.3333]);
ax = axes('parent', 2);
 

ax2 = axes('parent', 2, 'position', [0    0.9232    1    0.0623]);
t = text(ax2, .5, 0.5, 'Input 3D model', 'fontsize', 14, 'horizontalalignment', 'center');
set(ax2, 'visible', 'off');
axes(ax);
set(gca, 'color', 'none');
pobj = patch(obj,'FaceColor',       [0.6 .8 1.0], ...
    'EdgeColor',       'none',        ...
    'FaceLighting',    'gouraud',     ...
    'AmbientStrength', 0.15);
    
    % Add a camera light, and tone down the specular highlighting
camlight('headlight');
material('dull');

% Fix the axes scaling, and set a nice view angle
axis('image');
view([12 12 10]);
lightHand = findobj('parent', gca, 'type', 'light');

xlabel('Position (nm)', 'fontsize', 12);
ylabel('Position (nm)', 'fontsize', 12);
zlabel('Position (nm)', 'fontsize', 12);
%%
p =1;
for k = 120:3:359
    
    view(k+30, 19);
    axis([1000 1e4 1000 1e4 2000 1e4]);
    
    for m = 1:numel(lightHand)
        lightangle(lightHand(m), 60+k, 30)
    end
   
    drawnow;
    
          % Capture the plot as an image 
      frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if p == 1 
          imwrite(imind,cm,fileName,'gif', 'Loopcount',inf, 'DelayTime', gifDelay); 
          p = 0;
      else 
          imwrite(imind,cm,fileName,'gif','WriteMode','append', 'DelayTime', gifDelay); 
      end 
    
    
end

% Sampled point cloud
hold on
set(t, 'string', 'Generated SMLM point cloud');

ptsCld = plot3(pts(:,1), pts(:,2), pts(:,3), '.', 'markersize', 1, 'markeredgecolor', [0.8, 0.8, 0.8]);


for k = 1:3:179
    
    view(k+30, 19);
    axis([1000 1e4 1000 1e4 2000 1e4]);
    
	for m = 1:numel(lightHand)
        lightangle(lightHand(m), 60+k, 30)
    end
    
    
    drawnow;
    
    
              % Capture the plot as an image 
      frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 

          imwrite(imind,cm,fileName,'gif','WriteMode','append', 'DelayTime', gifDelay); 

    
end

set(pobj, 'visible', 'off');

for k = 180:3:359
    
    view(k+30, 19);
    axis([1000 1e4 1000 1e4 2000 1e4]);
    
    drawnow;
    
        
              % Capture the plot as an image 
      frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 

          imwrite(imind,cm,fileName,'gif','WriteMode','append', 'DelayTime', gifDelay); 

    
end

%%
% Add in fit curves
cla;

set(t, 'string', 'Piecewise point cloud fitting');

ptsCld = plot3(pts(:,1), pts(:,2), pts(:,3), '.', 'markersize', 1, 'markeredgecolor', [0.8, 0.8, 0.8]);

for k = 1:length(xPiecewise)
    
    if ~isempty(xPiecewise{k, 1})
    
        X = xPiecewise{k,2}*ones(size(xPiecewise{k}, 1), 1);
        Y = xPiecewise{k,1}(:,1);
        Z = xPiecewise{k,1}(:,2);
    
        Y = Y*cos(-rotationAroundX) - Z*sin(-rotationAroundX);
        Z = Y*sin(-rotationAroundX) + Z*cos(-rotationAroundX)+100;
    
        plot3(X, ...
            Y, ...
            1.1*Z, ...    
            'r', 'linewidth', 1);
    end
end

for k = 1:length(yPiecewise)
    
    if ~isempty(yPiecewise{k, 1})
        
        X = yPiecewise{k,1}(:,1);
        Y = yPiecewise{k,2}*ones(size(yPiecewise{k}, 1), 1);
        Z = yPiecewise{k,1}(:,2);
    
        Y = Y*cos(-rotationAroundX) - Z*sin(-rotationAroundX);
        Z = Y*sin(-rotationAroundX) + Z*cos(-rotationAroundX)+100;

        plot3(X, Y, 1.1*Z, ...    
            'b', 'linewidth', 1);
    end
end

for k = 1:3:179
    
    view(k+30, 19);
    axis([1000 1e4 1000 1e4 2000 1e4]);
    
    
    
    drawnow;
    
        
              % Capture the plot as an image 
      frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 

          imwrite(imind,cm,fileName,'gif','WriteMode','append', 'DelayTime', gifDelay); 

    
end

set(findobj('color', 'b'), 'visible', 'off');
set(findobj('color', 'r'), 'visible', 'off');

%%
% Add in fit surface

set(t, 'string', 'Noise-tolerant surface reconstruction');

tri = triangulation(plyFile.faces(:,2:4)+1, plyFile.vertices(:,1), plyFile.vertices(:,2), plyFile.vertices(:,3));
triPlot = trisurf(tri, 'edgecolor', 'none');
set(triPlot, 'faceColor', [0.6 1 0.4]);
material('dull');
camlight('headlight');
set(triPlot, ...
    'FaceLighting',    'gouraud', ...
    'AmbientStrength', 0.001)
lightHand = findobj('parent', gca, 'type', 'light');

set(ptsCld, 'visible', 'off');

for k = 180:3:269
    
    view(k+30, 19);
    axis([1000 1e4 1000 1e4 2000 1e4]);
    
	for m = 1:numel(lightHand)
        lightangle(lightHand(m), 60+k, 30)
    end
    
    
    drawnow;
    
        
              % Capture the plot as an image 
      frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 

          imwrite(imind,cm,fileName,'gif','WriteMode','append', 'DelayTime', gifDelay); 

    
end




% for k = 270:3:500
%     
%     view(k+30, 19);
%     axis([1000 1e4 1000 1e4 2000 1e4]);
%     
% 	for m = 1:numel(lightHand)
%         lightangle(lightHand(m), 60+k, 30)
%     end
%     
%     
%     drawnow;
%     
%         
%               % Capture the plot as an image 
%       frame = getframe(h); 
%       im = frame2im(frame); 
%       [imind,cm] = rgb2ind(im,256); 
%       % Write to the GIF File 
% 
%           imwrite(imind,cm,fileName,'gif','WriteMode','append', 'DelayTime', gifDelay); 
% 
%     
% end



set(t, 'string', 'Reconstructed surface curvature rendering');



set(gca,'CLim', [-.1 0.8]); % Set color limits of plot to 5th and 95th percentiles of angle measured
set(triPlot,'FaceColor','flat',...
       'FaceVertexCData',hsvColors(:,1),...
       'CDataMapping','scaled');
colormap('parula')

set(triPlot, ...
    'FaceLighting',    'gouraud', ...
    'AmbientStrength', 0.005)

for k = 270:3:469
    
    view(k+30, 19);
    axis([1000 1e4 1000 1e4 2000 1e4]);
    
	for m = 1:numel(lightHand)
        lightangle(lightHand(m), 60+k, 30)
    end
    
    
    drawnow;
    
        
              % Capture the plot as an image 
      frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 

          imwrite(imind,cm,fileName,'gif','WriteMode','append', 'DelayTime', gifDelay); 

    
end

% set(triPlot, 'visible', 'off');