function pixels = Lv(inpic, shape)
if (nargin < 2)
    shape = 'same';
end

dx = [0     -0.5    0
      0     0       0
      0     0.5     0];

dy = dx';

Lx = filter2(dx, inpic, shape);
Ly = filter2(dy, inpic, shape);
% pixels = Lx.^2 + Ly.^2;
pixels = sqrt(Lx.^2 + Ly.^2);
end