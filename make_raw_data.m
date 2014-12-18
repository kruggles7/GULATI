function [RBC_mean, RBC_std, parasite_mean, parasite_std, control] = make_raw_data( parasite, rbc, uninfected, names )

[r,c]=size(parasite); 

L=length(names); 
[r,c]=size(parasite);  
mat=parasite; 
mol_parasite=zeros(L,c); 
sum_mat=sum(parasite); 

%find the mol%time for each lipid group
for k=1:L;
    for j=1:c
        mol_parasite(k,j)=parasite(k,j)/sum_mat(1,j); 
    end 
end

[r,c]=size(rbc);  
mol_rbc=zeros(L,c); 
sum_mat=sum(rbc); 

%find the mol%time for each lipid group
for k=1:L
    for j=1:c
        mol_rbc(k,j)=rbc(k,j)/sum_mat(1,j);
    end 
end

[r,c]=size(uninfected);  
mol_uninfected=zeros(L,c); 
sum_mat=sum(uninfected); 

%find the mol%time for each lipid group
for k=1:L
    for j=1:c
        mol_uninfected(k,j)=uninfected(k,j)/sum_mat(1,j); 
    end 
end

RBC_mean=zeros(L,6); 
RBC_std=zeros(L,6); 
parasite_mean=zeros(L,6); 
parasite_std=zeros(L,6); 
control=zeros(L,2); 

for k=1:L
    c=1; 
    d=1;
    control(k,1)=mean(mol_uninfected(k,:))*100;
    std_dev=std(mol_uninfected(k,:)); 
    control(k,2)=(std_dev/sqrt(3))*100; 
    for j=1:6
        RBC_mean(k,j)=mean(mol_rbc(k,d:d+8))*100; 
        std_dev=std(mol_rbc(k,d:d+8)); 
        RBC_std(k,j)=(std_dev/sqrt(9))*100; 
        c=c+1;
        parasite_mean(k,j)=mean(mol_parasite(k,d:d+8))*100; 
        std_dev=std(mol_parasite(k,d:d+8)); 
        parasite_std(k,j)=(std_dev/sqrt(9))*100; 
        d=d+9; 
    end
end 

end

