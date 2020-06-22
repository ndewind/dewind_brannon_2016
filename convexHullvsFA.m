function [r] = convexHullvsFA

try
    tic
    % set parameters
    sess.pixelBuffer = 1;
    sess.numStim = 10000;
    
    % load stimuli
    load('559_stim.mat')
    
    % open onscreen window
    sca;
    [w sr] = Screen('OpenWindow',1,[0 0 0]);
    origin = [sr(3) sr(4)]/2;
    originRect = [origin origin];
    
    % preallocate
    stimIndx = nan(sess.numStim,1);
    CH = nan(sess.numStim,1);
    FA = nan(sess.numStim,1);
    
    for k = 1:sess.numStim
        % pick a random stimulus parameter set and instantiate it
        stimIndx(k) = ceil(rand*numel(allstim.tnum));
        pts = dotFieldGKA(round(allstim.num(stimIndx(k))),... % numerosity
            round(sqrt(allstim.fa(stimIndx(k))/pi)*2),... % field diameter
            round(sqrt(allstim.ia(stimIndx(k))/pi)*2) + sess.pixelBuffer); % distance btwn dot centers
        dots = pts2rect(pts,... % pts
            sqrt(allstim.ia(stimIndx(k))/pi),... % radius (can be vector if we do heterogenious dots)
            originRect,... % center of array in rect format (either left or right from coin flip above)
            round(sqrt(allstim.fa(stimIndx(k))/pi))); % radius of field
        
        % draw it and grab the image
        Screen('FillOval',w,[255 255 255],dots');
        Screen('Flip',w);
        dotArray = Screen('GetImage',w);
        
        % get coordinates of white pixels
        count = 0;
        imageCoord = nan(sum(dotArray(:))/255/3,2);
        for n = 1:size(dotArray,1)
            for m = 1:size(dotArray,2)
                if dotArray(n,m,1)
                    count = count+1;
                    imageCoord(count,1) = n;
                    imageCoord(count,2) = m;
                end
            end
        end
        
        % calculate convex hull
        [chIndx,CH(k)] = convhull(imageCoord);
        FA(k) = allstim.fa(stimIndx(k));
        
        % plot convex hull around dots (sanity check)
        % plot(imageCoord(chIndx,1),imageCoord(chIndx,2),'r-',imageCoord(:,1),imageCoord(:,2),'b.')
    end
    
    % cleanup
    sca
    [r,p]=corrcoef(CH,FA)
    toc
    keyboard
    
catch ME
    sca
    keyboard
end
end

function dots = pts2rect(pts,dotRads,centerRect,fieldRad)
if numel(dotRads) == 1
    dotRads = dotRads * ones(numel(pts)/2,1);
elseif numel(dotRads) ~= numel(pts)/2
    error('wrong number of dot rads')
end
fieldRadRect = [fieldRad fieldRad fieldRad fieldRad];
dots = [pts pts] + [-dotRads -dotRads dotRads dotRads] +...
    repmat(centerRect,numel(dotRads),1) -...
    repmat(fieldRadRect,numel(dotRads),1);
end
