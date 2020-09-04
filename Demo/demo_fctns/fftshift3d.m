% Fourier convention shift function
% (shift along 3 spatial dimensions, at each time)

function [A] = fftshift3d(A)

A = fftshift(fftshift(fftshift(A,1),2),3);

end