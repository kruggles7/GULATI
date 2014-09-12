function [ttest_mat, plot_mat] = Figure1_gams_no_infected( early, late, g_names, species )

cd ..
cd results
mkdir ('gams_bar_charts')
cd gams_bar_charts

g_groups=cell2mat(g_names(:,2)); 
unique_=unique(g_groups); 
indx=find(unique_~=0); 
unique_=unique_(indx); 
cmap=[255/255 51/255 51/255; 0/255 204/255 102/255; 204/255 0/255 204/255];
L_=length(unique_); 

[r,c]=size(early);  
mol_early=zeros(L_,c); 
sum_early=sum(early); 
 
mol_late=zeros(L_,c); 
sum_late=sum(late); 

%find the mol%time for each lipid group
for k=1:length(unique_);
    indx=find(g_groups==unique_(k)); 
    mat_e=early(indx,:); 
    mat_species_e=sum(mat_e); 
    mat_l=late(indx,:); 
    mat_species_l=sum(mat_l); 
   % for j=1:c
        mol_early(k,:)=mat_species_e; %/sum_early(1,j); 
        mol_late(k,:)=mat_species_l;%/sum_late(1,j); 
    %end 
    
%     indx=find(c1_groups==unique_(k)); 
%     mat_c1=control1(indx,:); 
%     mat_species_c1=sum(mat_c1); 
%     indx=find(c2_groups==unique_(k)); 
%     mat_c2=control2(indx,:); 
%     mat_species_c2=sum(mat_c2); 
%     for j=1:3
%         mol_t4(k,j)=mat_species_c1(1,j)/sum_c1(1,j); 
%         mol_t4(k,j+3)=mat_species_c2(1,j)/sum_c2(1,j); 
%     end 
end
ttest_mat=cell(length(unique_)+1,8); 

ttest_mat{1,2}='E vs. L parasite'; 
ttest_mat{1,3}='E vs. L iRBC'; 
ttest_mat{1,4}='E vs. L uiRBC'; 
ttest_mat{1,5}='E-gam vs. E-iRBC'; 
ttest_mat{1,6}='E-gam vs. E-uiRBC'; 
ttest_mat{1,7}='L-gam vs. L-iRBC'; 
ttest_mat{1,8}='L-gam vs. L-uiRNC'; 
ttest_mat{1,9}='E-iRBC vs. E-uiRBC'; 
ttest_mat{1,10}='L-iRBC vs. L-uiRBC'; 

plot_mat=zeros(L_,7); 
plot_std=NaN(L_,7); 
mol_names=cell(length(unique_),1); 

l_malaria=mol_late(:,1:9); 
l_infected=mol_late(:,10:18); 
l_uninf=mol_late(:,19:27); 
e_malaria=mol_early(:,1:9); 
e_infected=mol_early(:,10:18); 
e_uninf=mol_early(:,19:27); 

for k=1:length(unique_); 
    label=species{unique_(k),2}; 
    mol_names{k,1}=label; 
    ttest_mat{k+1,1}=label; 
    c=1; 
    d=1;
    for j=1:3
        if j==1 %gametocytes
            plot_mat(k,2)=mean(mol_early(k,d:d+8)); %*100; 
            std_dev=std(mol_early(k,d:d+8)); 
            plot_std(k,2)=(std_dev/sqrt(9)); %*100; 
            plot_mat(k,5)=mean(mol_late(k,d:d+8)); %*100; 
            std_dev=std(mol_late(k,d:d+8)); 
            plot_std(k,5)=(std_dev/sqrt(9)); %*100; 
        elseif j==2 %iRBC
            plot_mat(k,1)=mean(mol_early(k,d:d+8)); %*100; 
            std_dev=std(mol_early(k,d:d+8)); 
            plot_std(k,1)=(std_dev/sqrt(9)); %*100; 
            plot_mat(k,4)=mean(mol_late(k,d:d+8)); %*100; 
            std_dev=std(mol_late(k,d:d+8)); 
            plot_std(k,4)=(std_dev/sqrt(9)); %*100;   
        else 
            u=cat(2, mol_early(k,d:d+8), mol_late(k,d:d+8)) ;
            uiRBC(k,1)=mean(u) ;
            std_uiRBC=std(u);
            uiRBC(k,2)=std_uiRBC/sqrt(18); 
        end 
        [h,p]=ttest(mol_early(k,d:d+8),mol_late(k,d:d+8)); 
        ttest_mat{k+1,j+1}=p; 
        d=d+9; 
    end
    
    [h,p]=ttest(e_malaria(k,:), e_infected(k,:)); 
    ttest_mat{k+1,5}=p; 
    [h,p]=ttest(e_malaria(k,:), e_uninf(k,:)); 
    ttest_mat{k+1,6}=p; 
    [h,p]=ttest(l_malaria(k,:), l_infected(k,:)); 
    ttest_mat{k+1,7}=p; 
    [h,p]=ttest(l_malaria(k,:), l_uninf(k,:)); 
    ttest_mat{k+1,8}=p; 
    [h,p]=ttest(e_infected(k,:), e_uninf(k,:)); 
    ttest_mat{k+1,9}=p; 
    [h,p]=ttest(l_infected(k,:), l_uninf(k,:)); 
    ttest_mat{k+1,10}=p; 
    
    colors=[0 0 0; 0.5 0.5 0.5; 1 1 1]; 
    bar_color=[2 1 1 2 1]; 
    y=1:5; 
    for t=1:numel(y)
        bar( y(t), plot_mat(k,t),1,'facecolor',colors(bar_color(t),:)); 
        hold on 
        h(t)=errorbar (y(t),plot_mat(k,t),plot_std(k,t)); %,'.k', 'MarkerSize',2, 'linewidth',2 ); 
    end 
    
    % Perform the correction 
    for ii=1:length(h) 
            hc = (get(h(ii),'Children'))'; 
                % Get the x-coordinates for all children 
            xdata = get(hc(2),'XData'); 
                % Build arrays to select the proper x-coord. for all whiskers 
                % temp contains all indeces for the left x-coord of the whiskers 
            temp = 4:3:length(xdata); 
            temp(3:3:end) = []; 
                % temp2 contains the indeces for the original data point (where the whiskers are drawn around) 
            temp2=sort([1:9:length(xdata), 2:9:length(xdata)]); 
                % Build vectors containing indeces for the left and right x- coordinates of the whiskers 
            xleft = temp; xright = temp+1; 
         % Write the new data for the x-coordinates of the whiskers (original data point +- 0.01) 
            xdata(xleft) = xdata(temp2) - 0.01; 
            xdata(xright) = xdata(temp2) + 0.01; 
                % Set the new coordinates in the graph 
            set(hc(2),'Xdata',xdata, 'Color', 'k', 'LineWidth', 2) 
     end
    
    title(label, 'fontsize', 18); 
    x_control=0:10:30; 
    control=uiRBC(k,1); 
    control_up=uiRBC(k,1)+uiRBC(k,2);  
    control_down=uiRBC(k,1)-uiRBC(k,2); 
    y_control=[control control control control]; 
    y_control_up=[control_up control_up control_up control_up]; 
    y_control_down=[control_down control_down control_down control_down]; 
    plot(x_control, y_control, 'r-', 'linewidth', 2); 
    plot(x_control, y_control_up, 'r--', 'linewidth',2); 
    plot(x_control, y_control_down, 'r--', 'linewidth', 2); 
    
    
    hold off 
    %x_label= {'Gametocyte','iRBC', [],'Gametocyte','iRBC'}; 
    x_label={'Early', 'Late'}; 
    set(gca, 'XTick', [1.5 4.5]); %iRBC gray, gametocyte black
    set (gca, 'XTickLabel',x_label); 
    xlim([0 6]); 
    set(gca,'FontSize',18, 'linewidth',2);
    %NumTicks = 4;
    %L = get(gca,'YLim');
    %set(gca,'YTick',linspace(L(1),L(2),NumTicks))
    %ylabel('Mol % of Total Lipids Measured'); 
    print (gcf, '-depsc2', [label '_gametocyte_bar.eps']); 
    print (gcf, '-dpng', [label '_gametocyte_bar.png']); 
    close all 
end 

cd ..
cd ..
cd programs
end

