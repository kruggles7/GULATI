function [ pvalue_1_2, pvalue_1_3, pvalue_2_3 ] = ttest_within_species( data1, data2, data3, names, T)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%T=time
%S=samples

[r,c]=size(data1); 
%pvalue_mat = cell.empty;
pvalue_1_2(:,1)=names(:,1);
pvalue_1_3(:,1)=names(:,1);
pvalue_2_3(:,1)=names(:,1);

for i=1:r
        d=1; 
        time_3 = data3(i,:); 
        for j=1:T 
            time_1=data1(i,d:d+8); 
            time_2=data2(i,d:d+8); 
            [h,p]=ttest2(time_1, time_2); 
            pvalue_1_2 {i,j+1} = num2str(p); 
            [h,p]=ttest2(time_1, time_3,0.05,'both', 'unequal'); 
            pvalue_1_3{i,j+1}=num2str(p); 
            [h,p]=ttest2(time_2,time_3,0.05,'both', 'unequal');  
            pvalue_2_3{i,j+1}=p; 
            d=d+9; 
        end
    close all; 
end 

end

