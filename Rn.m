function V = Rn(n,gamma,nk)
    N= n;
    lammda = 3;
    t=length(nk);   
        % n! * k_t / k^n 
        v=0;
        for k = t:n
%             coeff=beta(k*gamma,N);
%             coeff = numer(t,k,coeff); 
%             v = v + coeff * poisspdf(k-1,lammda);
            coeff = betaln(k*gamma,N);
            coeff = coeff+gammaln(k+1)-gammaln(k-t+1);
            v = v + exp(coeff + (k-1)*log(lammda) - lammda - gammaln(k));
        end
    cnum=numel(nk);
    Be=0;
    for j = 1:cnum
        Be=Be-log(beta(gamma,nk(j)))+sum( log(   1:(nk(j)-1)   ) );
    end
    
V=log(v)+Be;

