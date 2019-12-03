function pixels = Lvvtilde(inpic, shape)
    if (nargin < 2)
        shape = 'same';
    end

    dx = [0     0       0       0       0;
          0     0       -0.5    0       0;
          0     0       0       0       0;
          0     0       0.5     0       0;
          0     0       0       0       0];
      
    dy = dx';
    
    dxx = [0     0       0       0       0;
           0     0       1       0       0;
           0     0       -2      0       0;
           0     0       1       0       0;
           0     0       0       0       0];
       
    dyy = dxx';
    
    dxy = conv2(dx, dy, 'same');
    
    Lx = conv2(inpic, dx, shape);
    Ly = conv2(inpic, dy, shape);
    
    Lxx = conv2(inpic, dxx, shape);
    Lxy = conv2(inpic, dxy, shape);
    Lyy = conv2(inpic, dyy, shape);
    
    pixels = (Lx.^2).*Lxx + 2.*Lx.*Ly.*Lxy + (Ly.^2).*Lyy;
end