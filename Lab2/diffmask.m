function [dx, dy] = diffmask(str)
    if strcmpi('central', str)
        dx = [0 -0.5 0; 0 0 0; 0 0.5 0];
        dx = [0     0   0;
              -0.5  0   0.5;
              0     0   0];
        dy = dx';
    elseif strcmpi('sobel', str)
        dx = [-1 -2 -1; 0 0 0; 1 2 1];
        dy = dx';
    end
end

