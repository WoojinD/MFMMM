function nprob = FH(X,c,n,k,gam,iii,kset,kset1,ksett,eachk,cij,ww) %c�� launch, each k �� k�� index

cc = unique(c); %component ���� ��
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


if isempty(find(ksett_k==iii))~=1 %claunch�� ci�� �󸶳� �ִ���, ����������� ����
%C_i �����ִ� X���� ���� ���ϱ�
i_=1;
for i = ksett% S U {i,j}
    if c(i)==cij %claunch�� �󺧰� ii�� �����͵��� ���� ����ִ´� 
        XC_igamma(:,i_)=Xgamma(:,i); %x����
        XC_igammalabel(:,i_)=i; %x���� ��ġ
        i_=i_+1;
    end
end
XC_igamma_k=XC_igamma;
k_label=find(eachk==XC_igammalabel);
if isempty(k_label)==0
XC_igamma_k(:,k_label)=[];
end

%mu0_gamma ���ϱ�
s=size(Xgamma);
ss=s(1);


%% kset ���� ��������
% j=1;
% for i = kset
% ksettt(j)=c(i);
% j=j+1;
% end


%%
%���⼭���� ����
XCS=size(XC_igamma);
XCS_k=size(XC_igamma_k);


%i,j���Ե� kset�� ���տ��� c(each)�� �����ֵ�  n-k,c is the number of cg for g =8 k in S U {i, j } that are eq
n_ci=XCS(2);
n_ci_i=XCS_k(2); 
L_gam=0;




    nprob = (n_ci_i)*exp(log(nprob)+detqs+detqs2);
else
    nprob=0;
end

end