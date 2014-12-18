function [ ttest_mat, N] = ttest_gametocyte_final( early, late, labels, cont_mat, cont_label)

[r,c]= size(early); 
ttest_mat = cell (r+1, 4); 

ttest_mat{1,2}='Early-control'; 
ttest_mat{1,3}='Late-control'; 
ttest_mat{1,4}='Early-Late'; 
N=nan(r,1); 

for j=1:r 
    name=labels{j};  
    E=early(j,:); 
    L=late(j,:); 
    T4=nan(9,1); 
    control=0; 
    indx=find(strcmp(cont_label,name)==1); 
    if numel(indx)>0
        T4(1:3,:)=cont_mat(indx,:); 
        control=1; 
    end  
    if control==1
        ttest_mat{j+1,1}=name; 
        [h,p]=ttest2(E,T4); 
        ttest_mat{j+1,2}=p; 
        [h,p]=ttest2(L,T4); 
        ttest_mat{j+1,3}=p; 
        [h,p]=ttest2(E,L); 
        ttest_mat{j+1,4}=p; 
    else 
        ttest_mat{j+1}=name; 
        ttest_mat{j+1,2}=NaN; 
        ttest_mat{j+1,3}=NaN; 
        [h,p]=ttest2(E,L); 
        ttest_mat{j+1,4}=p; 
    end 
end 

end

