function [ plot_mat,name_group ] = gametocyte_heatmap_step2( early, late, g_names, control1, c_names1, control2, c_names2)
% [r,c]=size(early); 
% early=e_malaria(14:r,:); 
% late=l_malaria(14:r,:); 
% g_names=gametocyte_text(14:r,:); 
% control1=T4_mat; 
% control2=T4_mat2; 
% c_names1=time4_cell; 
% c_names2=time4_2_cell; 

[r,c]=size(early); 
c_mat=zeros(r,6); %new control matrix 

%make new control matrix ------------------------
for j=1:r 
    name=g_names{j,1};  
    indx=find(strcmp(c_names1(:,1),name)==1); 
    if numel(indx)>0
        temp_mat(j,1:3)=control1(indx,:); 
    end 
    indx=find(strcmp(c_names2,name)==1); 
    if numel(indx)>0
        temp_mat(j,4:6)=control2(indx,:); 
    end 
end 

[c_mat,mat]=filter_std(temp_mat,6,1); 
[r,c]=size(early);  

mol_early=early; 
mol_late=late; 
mol_control=c_mat; 

%foreach time
plot_mat=double.empty; 
plot_total=zeros(1,24); 

%for each experiment
for k=1:r
    plot_total(1,1:6)= mol_control(k,1:6); 
    plot_total(1,7:15)=mol_early(k,1:9);
    plot_total(1,16:24)=mol_late(k,1:9);
    M=median(plot_total); 
    E=zeros(9,1); 
    L=zeros(9,1);
    C=zeros(6,1); 
    for j=1:9 
        E(j)=mol_early(k,j)/M; 
        L(j)=mol_late(k,j)/M; 
    end 
    for j=1:3
        C(j)=mol_control(k,j)/M; 
    end 
    plot_mat(k,1)=mean(C); 
    plot_mat(k,2)=mean(E);
    plot_mat(k,3)=mean(L); 
end 

plot_mat=log2(plot_mat); 

indx=find(isnan(plot_mat(:,1))==0 | isnan(plot_mat(:,2))==0 | isnan(plot_mat(:,3))==0); 
plot_mat2=plot_mat(indx,:); 
labels=g_names(indx,2); 

indx=find(isnan(plot_mat2(:,1))==1); 
indx2=find(isnan(plot_mat2(:,2))==1); 
indx3=find(isnan(plot_mat2(:,3))==1); 

plot_mat2(indx,1)=0; 
plot_mat2(indx2,2)=0; 
plot_mat2(indx3,3)=0;

plots=plot_mat2; 
name_group=labels; 
name_group=cell2mat(name_group); 
x={'T4', 'Early','Late'}; 
%plots=flipud(plot_mat2);  
%name_group=flipud(labels); 
cmap=redgreencmap(256); 
%cmap='jet'; 
H=HeatMap(plots,'ColumnLabels',x,'RowLabels',g_names(indx,1), 'DisplayRange', 1, 'Colormap', cmap, 'Symmetric','true'); 
plot(H);
print (gcf, '-dpng','gametocyte_heatmap.png'); 
print (gcf, '-depsc2', [ 'gametocyte_heatmap.eps']); 

%saveas(gcf, 'gametocyte_heatmap2.fig'); 
close all; 

end

