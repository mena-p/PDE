function r = outer_radius(y,zi)
%OUTER_RADIUS Summary of this function goes here
%   Detailed explanation goes here

r=(.102*(y/zi)^(1/3))*(1 -(.25*(y/zi)))*zi;

if r < 10
    r = 10;
end

end

