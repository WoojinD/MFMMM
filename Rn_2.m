function V = Rn_2(n,gamma)
    N= n;
    T = n;
    v = zeros(T,n);
    lammda = 3;


    for t = 1:T
        j=1;
        % n! * k_t / k^n 
        for k = t:n
            coeff = betaln(k*gamma,N)+700; %+700; *n!
            coeff = coeff+gammaln(k+1)-gammaln(k-t+1);
%           coeff = numer(t,k,coeff);
      %      coeff + (k-1)*log(lammda) - lammda - gammaln(k)
            v(t,j) = v(t,j) + coeff + (k-1)*log(lammda) - lammda - gammaln(k);% poisspdf(k-1,lammda);
            j=j+1;
        end
    end

    vv=exp(v);    
    for t = 1:(T-1)
        r(t) = vv(t+1)/ vv(t);
    end
    
    
%원랜 이게 맞음 근데 결과 거의 비슷
% for i = 1:21025
%     vvv(i)=sum(vv(i,:));
% end
% for t = 1:(T-1)
%     r(t) = vvv(t+1)/ vvv(t);
% end
%원랜 이게맞음
    
%     for t = 1:(T-1)
%         r(t) = exp(sum(v(t+1,:))-sum(v(t,:)));
%     end
% 
%     for t = 1:(T-1)
%          r(t) =exp(v(t+1)-v(t));
%      end

V=r;