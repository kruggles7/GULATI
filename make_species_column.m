[r,c]=size(DD2f_names); 

DD2f_names2=cell(r,3); 
for i=1:r
    s=DD2f_names{i,2}; 
    DD2f_names2(i,1:2)=DD2f_names(i,1:2); 
    species1=cell2mat(species(:,1)); 
    species2=cell2mat(species(:,3)); 
    indx=find(s==species1);
    DD2f_names2{i,3}=species2(indx(1)); 
end 
    