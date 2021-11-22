function nn = numer(t,k,coeff)
nn=coeff;
    for i = 1: t
        nn = nn*(k-i+1);
    end
end