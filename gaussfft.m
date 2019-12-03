function C = gaussfft(pic, t)

    [XSize, YSize] = size(pic);

    if mod(XSize, 2) == 0
        XEnd = XSize/2 - 1;
        XVek = -XSize/2:XEnd;
    else
        XEnd = XSize/2;
        XVek = -XSize/2 + 0.5:XEnd;
    end

    if mod(YSize, 2) == 0
        YEnd = YSize/2 - 1;
        YVek = -YSize/2:YEnd;
    else
        YEnd = YSize/2;
        YVek = -YSize/2 + 0.5:YEnd;
    end

    [X, Y] = meshgrid(XVek, YVek);

    filter = (1/(2*pi*t)*exp(-(X.^2+Y.^2)/(2*t)));
    filter = filter/sum(filter(:));
    
    VarEst = variance(filter);
%     figure
%     showgrey(filter)
    filterHat = abs(fft2(filter));               % Fourier transform of Gaussian filter

    Hhat = fft2(pic);                               % Fourier transform of picture

    Chat = filterHat.*Hhat;                      % Multiplication

    C = ifft2(Chat);                   % Invert the resulting Fourier transform
end




