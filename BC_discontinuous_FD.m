function flux = BC_discontinuous_FD(x,y)

    if abs(y) <= 30
        flux = 400;
    else
        flux = 0;
    end
    %flux = 400;
end

