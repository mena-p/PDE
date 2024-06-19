function flux = BC_discontinuous_FD(x,y)

    if abs(y) <= 5
        flux = -400;
    else
        flux = 0;
    end
end

