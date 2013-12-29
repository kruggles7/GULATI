function [mol_parasite, mol_rbc, mol_uninfected, plot_mat ] = Figure1( parasite, rbc, uninfected, names, species )

% parasite=DD2f_p_;    
% rbc=DD2f_R_;
% uninfected=uninf; 
% names=DD2f_names; 


groups=cell2mat(names(:,2)); 
unique_=unique(groups); 
cmap=[255/255 51/255 51/255; 0/255 204/255 102/255; 204/255 0/255 204/255];
L_=length(unique_); 
[r,c]=size(parasite);  
mat=parasite; 
mol_parasite=zeros(L_,c); 
sum_mat=sum(parasite); 

%find the mol%time for each lipid group
for k=1:length(unique_);
    indx=find(groups==unique_(k)); 
    mat_=mat(indx,:); 
    mat_species=sum(mat_); 
    for j=1:c
        mol_parasite(k,j)=mat_species(1,j)/sum_mat(1,j); 
    end 
end

[r,c]=size(rbc);  
mat=rbc; 
mol_rbc=zeros(L_,c); 
sum_mat=sum(rbc); 

%find the mol%time for each lipid group
for k=1:length(unique_);
    indx=find(groups==unique_(k)); 
    mat_=mat(indx,:); 
    mat_species=sum(mat_); 
    for j=1:c
        mol_rbc(k,j)=mat_species(1,j)/sum_mat(1,j); 
    end 
end

[r,c]=size(uninfected);  
mat=uninfected; 
mol_uninfected=zeros(L_,c); 
sum_mat=sum(uninfected); 

%find the mol%time for each lipid group
for k=1:length(unique_);
    indx=find(groups==unique_(k)); 
    mat_=mat(indx,:); 
    mat_species=sum(mat_); 
    for j=1:c
        mol_uninfected(k,j)=mat_species(1,j)/sum_mat(1,j); 
    end 
end
plot_mat=zeros(L_,24); 
plot_std=NaN(L_,24); 
for k=1:length(unique_); 
    label=species(k,2); 
    c=1; 
    d=1;
    for j=1:6
        plot_mat(k,c)=mean(mol_uninfected(k,:))*100;
        std_dev=std(mol_uninfected(k,:)); 
        plot_std(k,c)=(std_dev/sqrt(3))*100; 
        c=c+1;  
        plot_mat(k,c)=mean(mol_rbc(k,d:d+8))*100; 
        std_dev=std(mol_rbc(k,d:d+8)); 
        plot_std(k,c)=(std_dev/sqrt(9))*100; 
        c=c+1;
        plot_mat(k,c)=mean(mol_parasite(k,d:d+8))*100; 
        std_dev=std(mol_parasite(k,d:d+8)); 
        plot_std(k,c)=(std_dev/sqrt(9))*100; 
        c=c+2; 
        d=d+9; 
    end
    bar_color=[-1 0 1 0 -1 0 1 0 -1 0 1 0 -1 0 1 0 -1 0 1 0 -1 0 1 0]; 
    bh=bar(plot_mat(k,:), 1); 
    applyhatch(gcf,'\-x.');
    title(label); 
    hold on
    c=1; 
%    for j=1:6
        errorbar (plot_mat(k,:), plot_std(k,:), '.k', 'MarkerSize', 2);
%         c=c+4; 
%     end 
    hold off 
    k_str=num2str(k); 
    ch=get(bh,'Children');  
    set (ch,'CData',bar_color) ; 
    colormap (cmap); 
    %x_label={'C', 'R', 'P', [], 'C', 'R', 'P', [], 'C', 'R', 'P', [], 'C', 'R', 'P', [], 'C', 'R', 'P', [], 'C', 'R', 'P', []}; 
    x_label= {[],'8',[],[],[],'16',[],[],[],'24',[],[],[],'32',[],[],[],'40',[],[],[],'48',[],[]}; 
    set(gca, 'XTick', [1:24]); 
    set (gca, 'XTickLabel',x_label); 
    xlim([0 25]); 
    ylabel('Mol % of Total Lipids Measured'); 
    print (gcf, '-dpng', [k_str 'figure1_bar.png']); 
    close all 
end 

end

