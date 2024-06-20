function flux = BC_total_flux_FD(height,position)
%BC_NORMAL_FD Summary of this function goes here
%   Detailed explanation goes here

flux = 0.13/(5*sqrt(2*pi)) * exp(-0.5*(position./5).^2);

end

