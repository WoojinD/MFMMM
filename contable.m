function [conmatrix] = contable(label, result)

assert(length(label) == length(result));

label = label(:);
result = result(:);

n = length(label);

label_unique = unique(label);
result_unique = unique(result);

cls_num = max(label_unique);
res_num = max(result_unique); %numel(result_unique);

conmatrix = zeros(cls_num, res_num);

% compute the contigence table here

for i = 1: n
    rid = label(i);
    cid = result(i);
    newid = find(cid==result_unique);
    conmatrix(rid,cid) = conmatrix(rid,newid) + 1;
    %conmatrix(rid,cid) = conmatrix(rid,cid) + 1;
end

