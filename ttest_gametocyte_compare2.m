function [ ttest_mat, N] = ttest_gametocyte_compare2( early, late, labels)

[r,c]= size(early); 
ttest_mat = cell (r+1, 4); 

ttest_mat{1,2}='Early-T4'; 
ttest_mat{1,3}='Late-T4'; 
ttest_mat{1,4}='Early-Late'; 
N=nan(r,1); 

for j=1:r 
    name=labels{j};  
    E=early(j,:); 
    L=late(j,:); 
    ttest_mat{j+1}=name; 
    [h,p]=ttest2(E,L); 
    ttest_mat{j+1,2}=p; 
end 

end

