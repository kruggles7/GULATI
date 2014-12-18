function [ ttest_mat] = ttest_gametocyte2(mat1, mat2, mat3, labels)

[r,c]= size(mat1); 
ttest_mat = cell (r+1, 4); 
 
ttest_mat{1,2}='malaria-infected'; 
ttest_mat{1,3}='malaria-uninfected'; 
ttest_mat{1,4}='infected-uninfected'; 

for j=1:r 
    name=labels{j};  
    m1=mat1(j,:); 
    m2=mat2(j,:); 
    m3=mat3(j,:); 
    
    ttest_mat{j+1,1}=name; 
    [h,p]=ttest2(m1,m2); 
    ttest_mat{j+1,2}=p; 
    [h,p]=ttest2(m1,m3); 
    ttest_mat{j+1,3}=p; 
    [h,p]=ttest2(m2,m3); 
    ttest_mat{j+1,4}=p;  
end 

end

