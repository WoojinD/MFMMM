function logprob = GammagivenXc(X,c,n,k,gam,ww,C)
cc = unique(c); %component ���� ��
cn = unique(c);
p=length(gam);
pr=sum(gam);
Xgamma=X;
Xgammac=X;


if p~=pr %% ��� �������� Ȱ��Ǵ°� �ƴϸ�
%     for i = 1:p
%         if gam(p+1-i)==0
%             Xgamma(p+1-i,:)=[];
%         end
%     end
%     for i = 1:p
%         if gam(p+1-i)==1
%             Xgammac(p+1-i,:)=[];
%         end
%     end
    findgam=find(gam==0);
    Xgamma(findgam,:)=[];

    findgam2=find(gam==1);
    Xgammac(findgam2,:)=[];

    ik=1;
    for i = cc %component ���� xgamma ������
        kk=1;
        for j =1:n
            if i==c(j)
            XC_gamma{1,ik}{:,kk}=Xgamma(:,j);
            kk=kk+1;
            end
        end
        XC_gamma{2,ik}=i;
        ik=ik+1;
    end    

    ik=1;
    for i = cc %component ���� xgammac ������
        kk=1;
        for j =1:n
            if i==c(j)
            XC_gammac{1,ik}{:,kk}=Xgammac(:,j);
            kk=kk+1;
            end
        end
        XC_gammac{2,i}=i;
        ik=ik+1;
    end 

    aa=1;
    for az = cc %componment ���� length �����ֱ�
    cn(2,aa)=length(find(c==az));
    aa=aa+1;
    end
    clear aa;


%% ���⼭����
        ic=1;
        pxtheta=0;
        wgamma=ww;
        wgamma( :, find(gam==0) ) = [];
        for i = cc
            www= wgamma(i,:);
            ii=find(cell2mat(XC_gamma(2,:))==i);
            XCG=cell2mat(XC_gamma{1,ii});
            sXCG=size(XCG);
            pxtheta=pxtheta+C*sum(www*XCG)-0.5*sXCG(2)*sum(www.^2);
        end
        logprob=pxtheta;
        %% X^c�κ� likelihood ���ؾ���
%

%         wgammac=ww;
%         wgammac( :, find(gam==1) ) = [];
%         for i = cc
%             www= wgammac(i,:);
%             ii=find(cell2mat(XC_gammac(2,:))==i);
%             XCG=cell2mat(XC_gammac{1,ii});
%             sXCG=size(XCG);
%             pxtheta=pxtheta+C*sum(www*XCG)-0.5*sXCG(2)*sum(www.^2);
%         end
        

else %��� ������ Ȱ��Ǹ�
      logprob=-10^300;
end        
end