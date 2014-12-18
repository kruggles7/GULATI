function [RBC_mean, RBC_std, parasite_mean, parasite_std, control] = make_raw_data_gams( early, late, g_names )

L=length(g_names); 

[r,c]=size(early);  
mol_early=early; 
mol_late=late; 

parasite_mean=zeros(L,2); 
RBC_mean=zeros(L,2); 
parasite_std=zeros(L,2); 
RBC_std=zeros(L,2);
control=zeros(L,2); 

for k=1:L
    %parasite
    parasite_mean(k,1)=mean(mol_early(k,1:9)); %*100; 
    std_dev=std(mol_early(k,1:9)); 
    parasite_std(k,1)=(std_dev/sqrt(9)); %*100; 
    parasite_mean(k,2)=mean(mol_late(k,1:9)); %*100; 
    std_dev=std(mol_late(k,1:9)); 
    parasite_std(k,2)=(std_dev/sqrt(9)); %*100; 
    %iRBC
    RBC_mean(k,1)=mean(mol_early(k,10:18)); %*100; 
    std_dev=std(mol_early(k,10:18)); 
    RBC_std(k,1)=(std_dev/sqrt(9)); %*100; 
    RBC_mean(k,2)=mean(mol_late(k,10:18)); %*100; 
    std_dev=std(mol_late(k,10:18)); 
    RBC_std(k,2)=(std_dev/sqrt(9)); %*100; 
    %uiRBC
    u=cat(2, mol_early(k,19:27), mol_late(k,19:27)) ;
    control(k,1)=mean(u) ;
    std_uiRBC=std(u);
    control(k,2)=std_uiRBC/sqrt(18); 

end 

end

