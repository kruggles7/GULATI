function [mol_parasite, mol_rbc, mol_uninfected, plot_mat, plot_std, control_mat, control_std, mol_names, mol_parasite_all, mol_rbc_all, mol_uninf_all ] = Figure1_v3( parasite, rbc, uninfected, names, species )
cd ..
cd results
mkdir ('bar_charts')
cd bar_charts

[r,c]=size(parasite); 
mol_parasite_all=zeros(r,c); 
[r,c]=size(rbc); 
mol_rbc_all=zeros(r,c); 
[r,c]=size(uninfected); 
mol_uninf_all=zeros(r,c); 

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
        for i=1:numel(indx)
            mol_parasite_all(indx(i),j)=mat(indx(i),j)/sum_mat(1,j); 
        end 
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
        for i=1:numel(indx)
            mol_rbc_all(indx(i),j)=mat(indx(i),j)/sum_mat(1,j); 
        end 
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
        for i=1:numel(indx)
            mol_uninf_all(indx(i),j)=mat(indx(i),j)/sum_mat(1,j); 
        end 
    end 
end

plot_mat=zeros(L_,18); 
control_mat=zeros(L_,1); 
plot_std=NaN(L_,18); 
control_std=NaN(L_,1); 
mol_names=cell(length(unique_),1); 
for k=1:length(unique_); 
    label=species{unique_(k),2}; 
    mol_names{k,1}=label; 
    c=1; 
    d=1;
    control_mat(k,1)=mean(mol_uninfected(k,:))*100;
    std_dev=std(mol_uninfected(k,:)); 
    control_std(k,1)=(std_dev/sqrt(3))*100; 
    for j=1:6
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
    colors=[0 0 0; 0.5 0.5 0.5; 1 1 1]; 
    bar_color=[2 1 3 2 1 3 2 1 3 2 1 3 2 1 3 2 1 3]; 
    y=[1:18]; 
    for t=1:numel(y)
        bar( y(t), plot_mat(k,t),1,'facecolor',colors(bar_color(t),:)); %rbc=gray, parasite=black
        hold on 
        h(t)=errorbar (y(t),plot_mat(k,t),plot_std(k,t)); 
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
    control=control_mat(k,1); 
    control_up=control_mat(k,1)+control_std(k,1); 
    control_down=control_mat(k,1)-control_std(k,1); 
    y_control=[control control control control];
    y_control_up=[control_up control_up control_up control_up]; 
    y_control_down=[control_down control_down control_down control_down]; 
    
    plot(x_control, y_control, 'r-', 'linewidth', 2); 
    plot(x_control, y_control_up, 'r--', 'linewidth',2); 
    plot(x_control, y_control_down, 'r--', 'linewidth', 2); 
    
    hold off 
    x_label= {[],'8',[],[],'16',[],[],'24',[],[],'32',[],[],'40',[],[],'48',[]}; 
    set(gca, 'XTick', [1:18]); 
    set (gca, 'XTickLabel',x_label); 
    xlim([0 19]); 
    set(gca,'FontSize',18, 'linewidth',2);
    %NumTicks = 4;
    %L = get(gca,'YLim');
    %set(gca,'YTick',linspace(L(1),L(2),NumTicks))
    print (gcf, '-depsc2', [label '_figure1_bar.eps']); 
    print (gcf, '-dpng', [label '_figure1_bar.png']); 
    close all 
end 

cd ..
cd ..
cd programs

end

