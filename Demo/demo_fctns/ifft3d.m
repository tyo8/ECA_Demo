% Fourier transform function
% (transform along 3 spatial dimensions, at each time)

function [A] = ifft3d(A)

A = fftshift3d(ifft(ifft(ifft(ifftshift3d(A),[],1),[],2),[],3));

end