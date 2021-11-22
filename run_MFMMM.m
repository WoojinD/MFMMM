
load('erasedband.mat')
load('SalinasA_corrected.mat')
load('SalinasA.mat')
load('SalinasA_gt.mat')
clsid=reshape(salinasA_gt,[7138,1]);

data=reshape(salinasA,[7138,224]);
data=(data'-mean(mean(data)))'/mean(std(data));

for i = 1

iter=5000;

disp('dp_mmc clustering');
[params, tlast] = dpmm_mmc(data,iter, 4);
[mF] = F_measure(clsid(1:7138),params(end).classes);

lenpar=length(params);
for ii = 1:lenpar
[conmatrix] = contable(clsid(1:7138)'+1,params(ii).classes);
[v,hc,hk,h_ck,h_kc] = calculate_v_measure(conmatrix);
vhist(ii)=v;
end
[M,I]=max(vhist);
[mF_I] = F_measure(clsid(1:7138),params(I).classes);

%true image
figure
imshow(salinasA_gt,[],'InitialMagnification','fit');
title('true label image')

figure
data2=reshape(data,[83,86,224]);
data3=data2(:,:,10);
imshow(data3,[],'InitialMagnification','fit');


%clustering image
a=params(I).classes;
au=unique(a,'stable');
for i3 = 1:length(unique(a))
a(a==au(i3))=10+i3;
end
a=a-10;
aa=reshape(a,[83,86]);
figure
imshow(aa,[],'InitialMagnification','fit')
title('clustering image')


%band selected
selectedgam= find(params(iter+1).gam);
notselectedgam=setdiff([1:224],selectedgam);
difference=setdiff(erasedband,notselectedgam);
csvwrite('C:\Users\admin\Desktop\dpmmSVM\codeocean\difference.csv',difference) 

%score
csvwrite('C:\Users\admin\Desktop\dpmmSVM\codeocean\vscore.csv',M) 
csvwrite('C:\Users\admin\Desktop\dpmmSVM\codeocean\fscore.csv',mF_I) 
end


%%without backgroud
%%
idx=clsid;
for i = 1:7138
k=clsid(i);
if k ==0
idx(i)=0;
end
end

idxx=zeros(7138,1);
idxx(1)=1;
for i = 2:7138
    if ismember(idx(i),idxx)==1
        f=find(idx(i)==idxx);
        fi=f(1);
        idxx(i)=idxx(fi);
    else
        uni=unique(idxx);
        uninum=length(uni);
        idxx(i)=uninum;
    end
end

a=idxx;
aa=reshape(a,[83,86]);
figure
imshow(aa,[],'InitialMagnification','fit')
title('true label without background')

%%

idx=params(I).classes;
for i = 1:7138
k=clsid(i);
if k ==0
idx(i)=0;
end
end

idxx=zeros(7138,1);
idxx(1)=1;
for i = 2:7138
    if ismember(idx(i),idxx)==1
        f=find(idx(i)==idxx);
        fi=f(1);
        idxx(i)=idxx(fi);
    else
        uni=unique(idxx);
        uninum=length(uni);
        idxx(i)=uninum;
    end
end

a=idxx;
aa=reshape(a,[83,86]);
aa = uint8(255 * mat2gray(aa));
imwrite(aa,'C:\Users\admin\Desktop\dpmmSVM\codeocean\resultt.jpeg')
