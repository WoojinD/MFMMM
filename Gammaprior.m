function prob = Gammaprior(p,w,gamma)
prob = 0;
    for i = 1:p
    prob = prob+(gamma(i))*log(w)+(1-gamma(i))*log(1-w);
    end
end