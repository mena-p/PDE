function flux = BC_total_flux_FD(height,position)

flux = 0.13/(5*sqrt(2*pi)) * exp(-0.5*(position./5).^2);

end