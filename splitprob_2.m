function qL = splitprob_2(X,c,n,k,gam,iii,jjj,alpha,gamma,ccandi,ww,C)
%c(2)=2;%지울것
cc = unique(c); %component 종류 수
cn = unique(ccandi);
p=length(gam);
pr=sum(gam);
Xgamma=X;

wgamma=ww;
wgamma( :, find(gam==0) ) = [];

s=size(wgamma);
ws=s(1);


findgam=find(gam==0);
Xgamma(findgam,:)=[];
    
s=size(Xgamma);
ss=s(1);


%%p(C)
xs=size(X);
n=xs(2);
j=1;
for i = cc
ff=find(c==i);
nk(j)=length(ff);
j=j+1;
end

j=1;
for i = cn
ff=find(ccandi==i);
nk2(j)=length(ff);
j=j+1;
end


pC = Rn(n,gamma,nk);
pC_split = Rn(n,gamma,nk2);
%%
%여기서부터 시작
    Xiii=Xgamma(:,iii);
    ciii=c(iii);
    wwi= wgamma(ciii,:);
    
    Xjjj=Xgamma(:,jjj);
    cjjj=c(jjj);
    wwj= wgamma(cjjj,:);
    
    FGiii=C*wwi*Xiii-0.5*sum(wwi.^2);
    FGjjj=C*wwj*Xjjj-0.5*sum(wwj.^2);
    
    L_split= FGiii-FGjjj ;
  %  L_split = L_split+
    %분모-> 분자왼쪽-> 분자 오른쪽 순서
    %qL=exp(pC_split-pC)*exp(L_split)/(det( eye(pr)*k1+1/(1+h1)* (Xgamma(:,iii)-mu0_gamma)*(Xgamma(:,iii)-mu0_gamma)'  )*det( eye(pr)*k1+1/(1+h1)* (Xgamma(:,jjj)-mu0_gamma)*(Xgamma(:,jjj)-mu0_gamma)'))^((delta+pr)/2);
    
    qL=pC-pC_split+L_split;
    
    qL=exp(qL);
end