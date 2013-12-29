function [ DD2_new ] = reorder_DD2( DD2 )

[r,c]=size(DD2_mat); 

i=1; 
p=1; 
for k=1:6 %1,2,3,19,20,21,37,38,39
    for j=1:3
        DD2_new(:,p:p+2)=DD2(:,i:i+2); 
        p=p+18; 
        i=i+3; 
    end 
end 
end

