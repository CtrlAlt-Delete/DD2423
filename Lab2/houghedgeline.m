function [linepar acc] = houghedgeline(pic, scale, gradmagnthreshold, nrho, ntheta, nlines, verbose)
curves = extractedge(pic, scale, gradmagnthreshold, 'same');
smoothimg = discgaussfft(pic, scale);
magnitude = (Lv(smoothimg, 'same'));   % kanske ha sqrt över hela

[linepar, acc] = houghline(curves, magnitude, nrho, ntheta, gradmagnthreshold, nlines, verbose, pic, scale);
end