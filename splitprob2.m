function qL = splitprob2(X,c,gam,iii,jjj,alpha,prodfqb,csplitmerge,gamma,ww,C)
%c(2)=2;%지울것
cc = unique(c); %component 종류 수
cn = unique(csplitmerge);
p=length(gam);
pr=sum(gam);
Xgamma=X;

wgamma=ww;
wgamma( :, find(gam==0) ) = [];

findgam=find(gam==0);
Xgamma(findgam,:)=[];


s=size(Xgamma);
ss=s(1);
n=s(2);


%mu0_gamma 구하기
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
ff=find(csplitmerge==i);
nk2(j)=length(ff);
j=j+1;
end


pC = Rn(n,gamma,nk);
pC_split = Rn(n,gamma,nk2);


%%
%여기서부터 시작
    qccqcc=prodfqb; %%%%%%%%%나중에 다 exp 해야한다
    pcpc=pC_split-pC;

    FGiii=0;
    for ii=find(csplitmerge(iii)==csplitmerge)% k:csplit_k=csplit_i
        iiiwidx=csplitmerge(iii);
        try
        iiiw=wgamma(iiiwidx,:);
        catch
        iiiw=wgamma(end,:);
        end
        FGiii=FGiii+C*iiiw*Xgamma(:,ii)-0.5*sum(iiiw.^2);
    end
    
    FGjjj=0;
    for ii=find(csplitmerge(jjj)==csplitmerge)
        jjjwidx=csplitmerge(jjj);
        jjjw=wgamma(jjjwidx,:);
        FGjjj=FGjjj+C*jjjw*Xgamma(:,ii)-0.5*sum(jjjw.^2);
    end 
      
    FGccc=0;
    for ii=find(c(iii)==c)
        cccwidx=c(iii);
        cccw=wgamma(cccwidx,:);
        FGccc=FGccc+C*cccw*Xgamma(:,ii)-0.5*sum(cccw.^2);
    end    


    
    lik=FGiii+FGjjj-FGccc;
    qL=exp(qccqcc+pcpc+lik);

end