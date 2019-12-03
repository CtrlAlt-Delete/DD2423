function [linepar, acc] = houghline(curves, magnitude, nrho, ntheta, threshold, nlines, verbose, pic, scale)
% Check if input appear to be valid
% Allocate accumulator space
acc = zeros(nrho, ntheta);

% Define a coordinate system in the accumulator space
thetaRange = linspace(-pi/2, pi/2,ntheta);
[xD, yD] = size(magnitude);
D = sqrt(xD^2 + yD^2);
rhoRange = linspace(-D, D, nrho);

% Loop over all the input curves (cf. pixelplotcurves)
insize = size(curves, 2);   % Amount of edge pixels
trypointer = 1;
numcurves = 0;

% For each point on each curve
while trypointer <= insize
    polylength = curves(2, trypointer);
    numcurves = numcurves + 1;
    trypointer = trypointer + 1;
    
    for polyidx = 1:polylength
        x = curves(2, trypointer);
        y = curves(1, trypointer);
        
        % Check if valid point with respect to threshold
        if (magnitude(round(x), round(y)) > threshold)
            
            % Optionally, keep value from magnitude image
            magnValue = magnitude(round(x), round(y));
            
            % Loop over a set of theta values
            for thetaIndex = 1:length(thetaRange)
                
                % Compute rho for each theta value
                rhoVec(thetaIndex) = x*cos(thetaRange(thetaIndex)) + y*sin(thetaRange(thetaIndex));
                
                % Compute index values in the accumulator space
                [dist, rhoIndex] = min(abs(rhoRange-rhoVec(thetaIndex)));
                
                % Update the accumulator
%                 acc(rhoIndex, thetaIndex) = acc(rhoIndex, thetaIndex) + 1;            % use this line for standard voting
                acc(rhoIndex, thetaIndex) = acc(rhoIndex, thetaIndex) + log(magnValue); % use this line for weighted voting
            end
        end
        trypointer = trypointer + 1;
    end
end

% Extract local maxima from the accumulator
[pos value] = locmax8(acc);
[dummy indexvector] = sort(value);
nmaxima = size(value, 1);

% Delimit the number of responses if necessary
% Compute a line for each one of the strongest responses in the accumulator
linepar = [];
outcurves = zeros(2, 4*nlines);
for idx = 1:nlines
    rhoidxacc = pos(indexvector(nmaxima - idx + 1), 1);
    thetaidxacc = pos(indexvector(nmaxima - idx + 1), 2);
    rhoMax = rhoRange(rhoidxacc);
    thetaMax = thetaRange(thetaidxacc);
    linepar(:, idx) = [rhoMax; thetaMax];

    x0 = rhoMax*cos(thetaMax);       % some scaling issue? minus D? Center image?
    y0 = rhoMax*sin(thetaMax);
    dx = D^2*(-sin(thetaMax));
    dy = D^2*(cos(thetaMax));    
    
    outcurves(1, 4 * (idx - 1) + 1) = 0;
    outcurves(2, 4 * (idx - 1) + 1) = 3;
    outcurves(2, 4 * (idx - 1) + 2) = x0 - dx;
    outcurves(1, 4 * (idx - 1) + 2) = y0 - dy;
    outcurves(2, 4 * (idx - 1) + 3) = x0;
    outcurves(1, 4 * (idx - 1) + 3) = y0;
    outcurves(2, 4 * (idx - 1) + 4) = x0 + dx;
    outcurves(1, 4 * (idx - 1) + 4) = y0 + dy;
end

% Overlay these curves on the gradient magnitude image
figure
% subplot(1,2,1)
% overlaycurves(magnitude, outcurves);  % This line for testimages
overlaycurves(pic, outcurves);          % This line for real images
axis([1 size(magnitude, 2) 1 size(magnitude, 1)]);
% title('overlaycurves')    % This line for testimages
title(['scale = ', num2str(scale), ', threshold = ', num2str(threshold),...  % This line for real images
    ', nrho = ', num2str(nrho), ', ntheta = ', num2str(ntheta),...
    ', nlines = ', num2str(nlines)])
axis on
xlabel('y')
ylabel('x')

% Return the output data
end