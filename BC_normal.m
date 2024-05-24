function flux = BC(x)
%BC Summary of this function goes here
%   Detailed explanation goes here
sigma = 5;
flux = 0.01 * 1/(sigma*sqrt(2*pi)) * exp(-0.5.*(x./sigma).^2); % mol per second per sqr meter

end

