function T=tense_prod_index(v1,v2)

if isempty(v1)
    T=v2;
    return
end

if isempty(v2)
    T=v1;
    return
end


% v1 is a matrix and v2 is also matrix
% tensor prod along the index

[v1r,v1c]=size(v1);
[v2r,v2c]=size(v2);


T=padarray(repmat(v1,v2r,1),[0,v2c],0,'post');

for i=1:1:v2r
   T((i-1)*v1r+1:i*v1r,v1c+1:end)=repmat(v2(i,:),v1r,1); 
end




