function edgecurves = extractedge(inpic, scale, threshold, shape)

    smoothimg = discgaussfft(inpic, scale);
%     gradmagnimg = sqrt(Lv(smoothimg, 'same'));
    gradmagnimg = Lv(smoothimg, 'same'); 
    
    Lvv = Lvvtilde(discgaussfft(smoothimg, scale), shape);
    Lvvv = Lvvvtilde(discgaussfft(smoothimg, scale), shape);
    
    Lvvvmask = (Lvvv < 0) - 0.5;
    Lvmask = (gradmagnimg > threshold) - 0.5;
    
    edgecurves = zerocrosscurves(Lvv, Lvvvmask);
    edgecurves = thresholdcurves(edgecurves, Lvmask);
end