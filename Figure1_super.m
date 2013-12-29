function [mol_super,mol_control, mol_uninf, plot_mat, mol_names] = Figure1_super( super, control, names, species )

% parasite=DD2f_p_;    
% rbc=DD2f_R_;
% uninfected=uninf; 
% names=DD2f_names; 
% super=super_; 
% control=control_mean; 
% uninf=uiRBC_mean; 
% names=super_names; 

 %get information for controls 
[r,c]=size(control); 
control_new=zeros(r,9);  
control_mean=zeros(r,1); 
control_std=zeros(r,1); 
uiRBC_new=zeros(r,9); 
uiRBC_mean=zeros(r,1); 
uiRBC_std=zeros(r,1); 
for i=1:r
    %go through each of the time points
    timed=control(i,1:9); 
    means=mean(timed); 
    stds=std(timed); 
    filter_max=means+2*stds; 
    filter_min=means-2*stds; 
    indx=find(timed>filter_max | timed<filter_min); 
    timed(indx)=NaN; 
    control_new(i,:)=timed; 
    control_mean(i,1)=nanmean(timed); 
    n=isnan(timed); %how many are NaN;  
    if n<9
        control_std(i,1)=nanstd(timed)/sqrt(9-sum(n)); 
    else 
        control_std(i,1)=NaN; 
    end 
    
    timed_R=control(i,10:18);
    means=mean(timed_R); 
    stds=std(timed_R); 
    filter_max=means+2*stds; 
    filter_min=means-2*stds; 
    indx=find(timed_R>filter_max | timed_R<filter_min); 
    timed_R(indx)=NaN; 
    uiRBC_new(i,:)=timed_R; 
    uiRBC_mean(i,1)=nanmean(timed_R); 
    n=isnan(timed_R); %how many are NaN;  
    if n<9
        uiRBC_std(i,1)=nanstd(timed_R)/sqrt(9-sum(n)); 
    else 
        uiRBC_std(i,1)=NaN; 
    end 
    
end 

groups=cell2mat(names(:,2)); 
unique_=unique(groups); 
cmap=[255/255 51/255 51/255; 0/255 204/255 102/255; 204/255 0/255 204/255];
L_=length(unique_); 
[r,c]=size(super);  
mat=super; 
mol_super=zeros(L_,c); 
sum_mat=sum(super); 

%find the mol%time for each lipid group
for k=1:length(unique_);
    indx=find(groups==unique_(k));
    if numel(indx)>1
        mat_=mat(indx,:); 
        mat_species=sum(mat_); 
        for j=1:c
            mol_super(k,j)=mat_species(1,j)/sum_mat(1,j); 
        end 
    end 
end


[r,c]=size(control_new);  
mat=control_new; 
mol_control=zeros(L_,c); 
sum_mat=nansum(control_new); 

%find the mol%time for each lipid group
for k=1:length(unique_);
    indx=find(groups==unique_(k)); 
    if numel(indx)>1
        mat_=mat(indx,:); 
        mat_species=nansum(mat_);
        for j=1:c
            mol_control(k,j)=mat_species(1,j)/sum_mat(1,j); 
        end
    end 
end


[r,c]=size(uiRBC_new);  
mat=uiRBC_new; 
mol_uninf=zeros(L_,c); 
sum_mat=nansum(uiRBC_new); 

%find the mol%time for each lipid group
for k=1:length(unique_);
    indx=find(groups==unique_(k)); 
    if numel(indx)>1
        mat_=mat(indx,:); 
        mat_species=nansum(mat_); 
        for j=1:c
            mol_uninf(k,j)=mat_species(1,j)/sum_mat(1,j); 
        end
    end 
end

plot_mat=zeros(L_,6); 
plot_std=NaN(L_,6); 
mol_names=cell(length(unique_),1); 
for k=1:length(unique_); 
    label=species{unique_(k),2}; 
    mol_names{k,1}=label; 
    c=1; 
    d=1;
    C=nanmean(mol_control(k,:))*100;
    U=nanmean(mol_uninf(k,:))*100; 
    for j=1:6
        plot_mat(k,c)=mean(mol_super(k,d:d+8))*100; 
        std_dev=std(mol_super(k,d:d+8)); 
        plot_std(k,c)=(std_dev/sqrt(9))*100; 
        c=c+1; 
        d=d+9; 
    end
    bar_color=[-1 0 0 -1 0 0 -1 0 0 -1 0 0 -1 0 0 -1 0 0]; 
    bh=bar(plot_mat(k,:), 0.5, 'k'); 
    title(label); 
    hold on
    c=1; 
    errorbar (plot_mat(k,:), plot_std(k,:), '.k', 'MarkerSize', 2);
    x_control=0:10:30; 
    y_control=[C C C C]; 
    plot(x_control, y_control, 'r-'); 
    y_uninf=[U U U U];  
    plot(x_control, y_uninf, 'b--');  
    hold off 
    k_str=num2str(k); 
    x_label= {'8' '16' '24' '32' '40' '48'}; 
    set(gca, 'XTick', [1:6]); 
    set (gca, 'XTickLabel',x_label); 
    xlim([0 7]); 
    ylabel('Mol % of Total Lipids Measured'); 
    print (gcf, '-dpng', [label '_figure1_super.png']); 
    close all 
end 
%end

