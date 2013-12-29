function [ mat_filter, mat_NaN ] = filter_std( mat,s,t )
%filters every time point by average +- std 
%s=the number of samples in each timepoint 
[r,c]=size(mat); 
mean_=zeros(r,6); 
std_=zeros(r,6);  
c=1; 
mat_filter=zeros(r,c); 
mat_NaN=zeros(r,c);

for i=1:t
    for j=1:r
        e=c+s-1; 
        d=mat(j,c:e) ;
        m=mean(d);
        st=std(d); 
        mean_(j,i)=m; 
        std_(j,i)=st; 
    end 
    c=c+s; 
end 

c=1; 
for i=1:t
    for j=1:r
        filter_up=mean_(j,i)+2*std_(j,i); 
        filter_dn=mean_(j,i)-2*std_(j,i); 
        e=c+s-1;
        d=mat(j,c:e);
        d2=mat(j,c:e); 
        fup=filter_up; 
        fdn=filter_dn; 
        indx=find(d>fup);   
        d(indx)=NaN; 
        d2(indx)=mean_(j,i); 
        indx=find(d<fdn);   
        d(indx)=NaN; 
        d2(indx)=mean_(j,i); 
        mat_NaN(j,c:e)=d;
        mat_filter(j,c:e)=d2; 
    end 
    c=c+s; 
end 

end

