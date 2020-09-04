% Fourier convention shift function
% (shift along 3 spatial dimensions, at each time)

function [A] = ifftshift3d(A)

A = ifftshift(ifftshift(ifftshift(A,1),2),3);

end