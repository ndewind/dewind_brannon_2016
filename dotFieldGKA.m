function pts = dotFieldGKA(ndots, fieldSize, interDot)
% dotFieldGKA   Generate a random dot field
% Usage:
%     pts = dotFieldGKA(ndots, fieldSize, interDot)
%   Specify a number of points within a circular field that are positioned
% randomly, subject only to the constraint that they must be some distance
% from each other.  These points can be used as the centers of circular
% dots for a field of random, non-overlapping dots.
%   ndots is the number of points to generate.
%   fieldSize is the diameter of the field within which the points are to
% be generated, specified in pixels.
%   interDot is the minimum distance between points, specified in pixels.
% To produce non-overlapping dots, this value could be the diameter of the
% dots plus a small buffer.
%   pts is an ndots-by-2 matrix, such that pts(k,1) is the x coordinate of
% the k'th dot, and pts(k,2) is its y coordinate. Note that this is
% coordinates within the dot field; to convert to screen coordinates for
% drawing, the upper left-hand coordinate of the dot field should be added
% to the pts values.
%
% 2009/12/25 GKA
%
% interdot can now be a vector of length ndots (or still a simple scalar),
% which allows dots of different sizes within a single field. When a vector
% is used each value is the dot radius.  when a scalar is used interDot is
% still the diameter plus a buffer.

% Also replaced 'dummy' variable with a '~'

% 2010/7/27 NKD

if length(interDot) > 1 && length(interDot) ~= ndots
    error('interDot~=ndots','interdot is not scalar and does not equal ndots')
end

field = rand(fieldSize);
[Y, X] = meshgrid(0:(fieldSize-1));
rad = (fieldSize-1) / 2;
outField = sqrt((X-rad).^2 + (Y-rad).^2) > rad;
field(outField) = 0;

field = field(:);
X = X(:);
Y = Y(:);

pts = zeros(ndots, 2);

siz = [fieldSize fieldSize];

buffer = 2;

for k=1:ndots
    currField=field;
    if length(interDot) > 1
        for i=1:k-1
            tooClose = sqrt((X - (pts(i,1)-1)).^2 + (Y - (pts(i,2)-1)).^2)...
                <= interDot(k)+interDot(i)+buffer;
            currField(tooClose)=0;
        end
        [asdf, I] = max(currField);
        [xk, yk] = ind2sub(siz, I);
        pts(k,:) = [xk yk];
        if ~any(currField)
            error('dotFieldGKA:noSpace', 'Can''t fit in any more dots!');
            break;
        end
    else
        [asdf, I] = max(field);
        [xk, yk] = ind2sub(siz, I);
        tooClose = sqrt((X - (xk-1)).^2 + (Y - (yk-1)).^2) <= interDot;
        field(tooClose) = 0;
        pts(k,:) = [xk yk];
        if ~any(field)
            error('dotFieldGKA:noSpace', 'Can''t fit in any more dots!');
            break;
        end
    end
end
end