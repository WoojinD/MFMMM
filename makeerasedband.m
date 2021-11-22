load('SalinasA.mat')
load('SalinasA_corrected.mat')



k=1;
for i = 1: 220
    a =salinasA(:,:,i);
    c=0;
    for j = 1:200
        b=salinasA_corrected(:,:,j);
        if 1==isequal(a,b)
            c=1;
        end
    end
    if c==0
        erasedband(k)=i;
        k=k+1;
    end
end