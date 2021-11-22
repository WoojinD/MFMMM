function nprob = FH(X,c,n,k,gam,iii,kset,kset1,ksett,eachk,cij,ww) %c는 launch, each k 는 k의 index

cc = unique(c); %component 종류 수
cn = unique(c);
aa=1;
p=length(gam);
pr=sum(gam);
Xgamma=X;
Xgammac=X;

kkkkk=find(ksett==eachk);
ksett_k=ksett;
ksett_k(kkkkk)=[];

findgam=find(gam==0);
Xgamma(findgam,:)=[];


if isempty(find(ksett_k==iii))~=1 %claunch에 ci가 얼마나 있는지, 비어있을수도 있음
%C_i 갖고있는 X들의 집합 구하기
i_=1;
for i = ksett% S U {i,j}
    if c(i)==cij %claunch의 라벨과 ii가 같은것들은 여기 집어넣는다 
        XC_igamma(:,i_)=Xgamma(:,i); %x값들
        XC_igammalabel(:,i_)=i; %x값들 위치
        i_=i_+1;
    end
end
XC_igamma_k=XC_igamma;
k_label=find(eachk==XC_igammalabel);
if isempty(k_label)==0
XC_igamma_k(:,k_label)=[];
end

%mu0_gamma 구하기
s=size(Xgamma);
ss=s(1);


%% kset 집합 만들어야함
% j=1;
% for i = kset
% ksettt(j)=c(i);
% j=j+1;
% end


%%
%여기서부터 시작
XCS=size(XC_igamma);
XCS_k=size(XC_igamma_k);


%i,j포함된 kset에 집합에서 c(each)랑 같은애들  n-k,c is the number of cg for g =8 k in S U {i, j } that are eq
n_ci=XCS(2);
n_ci_i=XCS_k(2); 
L_gam=0;




    nprob = (n_ci_i)*exp(log(nprob)+detqs+detqs2);
else
    nprob=0;
end

end