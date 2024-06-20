function flux = BC_normal_FD(height,position)
%BC_NORMAL_FD Summary of this function goes here
%   Detailed explanation goes here

flux = -5175/(5*sqrt(2*pi)) * exp(-0.5*(position./5).^2);

end