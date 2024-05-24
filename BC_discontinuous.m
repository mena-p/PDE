function flux = BC_discontinuous(x)

flux = zeros(1,size(x,2));

for i = 1:size(x,2)
    if abs(x(1,i)) <= 5
        flux(i) = 0.01;
    else
        flux(i) = 0;
    end
end

end

