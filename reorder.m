function [ mat_new ] = reorder( mat )

[r,c]=size(mat); 
num=c/6;
n=sqrt(num); 
mat_new=zeros(r,c); 
i=1; 
p=1; 
r=1; 
for k=1:6 %1,2,3,19,20,21,37,38,39,4,5,6
    p=r; 
    for j=1:n
        mat_new(:,i:i+2)=mat(:,p:p+2); 
        p=p+18; 
        i=i+3; 
    end 
    r=r+3; 
end 
end

