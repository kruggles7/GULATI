function [ A ] = gametocyte_heatmap( early, late, g_names, control1, control2, c_names1, c_names2, species)

[r,c]=size(early); 
early=e_malaria(14:r,:); 
late=l_malaria(14:r,:); 
g_names=gametocyte_text(14:r,:); 
control1=time4; 
control2=time4_2; 
c_names1=time4_cell; 
c_names2=time4_2_cell; 

[r,c]=size(early); 
c_mat=zeros(r,6); %new control matrix 
c_mat2=zeros(r,3); 
temp_mat=zeros(r,6); 
temp_mat2=zeros(r,3); 

%make new control matrix ------------------------
for j=1:r 
    name=g_names{j,1};  
    indx=find(strcmp(c_names1(:,1),name)==1); 
    if numel(indx)>0
        temp_mat(j,1:3)=control1(indx,:); 
        temp_mat2(j,1:3)=control1(indx,:); 
    end 
    indx=find(strcmp(c_names2,name)==1); 
    if numel(indx)>0
        temp_mat(j,4:6)=control2(indx,:); 
    end 
end 

[c_mat,mat]=filter_std(temp_mat,6,1); 
%[c_mat,mat]=filter_std(temp_mat2,3,1); 
[r,c]=size(early);  
sum_early=sum(early); 
sum_late=sum(late); 
sum_c=sum(c_mat); 

mol_early=early; 
mol_late=late; 
mol_control=c_mat; 

%foreach time
plot_mat=double.empty; 
plot_total=zeros(1,24); 
%plot_total=zeros(1,21); 

%for each experiment
for k=1:r
    plot_total(1,1:6)=mol_control(k,1:6); 
    plot_total(1,7:15)=mol_early(k,1:9);
    plot_total(1,16:24)=mol_late(k,1:9); 
    %plot_total(1,1:3)=mol_control(k,1:3); 
    %plot_total(1,4:12)=mol_early(k,1:9); 
    %plot_total(1,13:21)=mol_early(k,1:9); 
    M=median(plot_total); 
    E=zeros(9,1); 
    L=zeros(9,1);
    %C=zeros(6,1); 
    C=zeros(3,1); 
    for j=1:9 
        E(j)=mol_early(k,j)/M; 
        L(j)=mol_late(k,j)/M; 
    end 
    for j=1:6
        C(j)=mol_control(k,j)/M; 
    end 
    plot_mat(k,1)=mean(C); 
    plot_mat(k,2)=mean(E);
    plot_mat(k,3)=mean(L); 
end
raw=plot_mat; 

indx=find(plot_mat==0); 
plot_mat(indx)=1; 
indx=find(isnan(plot_mat)==1); 
plot_mat(indx)=1; 
plot_mat=log2(plot_mat); 
x={'8','16','24','32','40','48'}; 
plots=flipud(plot_mat); 
%G=flipud(G); 
%name_group=names(:,1); 
cmap=redgreencmap(256); 
%cmap='jet'; 
x={'T4', 'early', 'late'}; 
%
H=HeatMap(plots,'ColumnLabels',x,'DisplayRange', 2, 'Colormap', cmap, 'Symmetric','true'); 
set(gca,'FontSize',16,'linewidth',2)
plot(H);
print (gcf, '-dpng', 'gametocyte_heatmap.png'); 
print (gcf, '-depsc2',  'gametocyte_heatmap.eps'); 
saveas(gcf, 'gametocyte_heatmap.fig'); 
close 


indx=find(isnan(plot_mat(:,1))==0 | isnan(plot_mat(:,2))==0 | isnan(plot_mat(:,3))==0); 
plot_mat2=plot_mat(indx,:); 
raw=raw(indx,:); 
CONTROL=mol_control(indx,:); 
EARLY=mol_early(indx,:); 
LATE=mol_late(indx,:); 
labels=g_names(indx,1); 

indx=find(isnan(plot_mat2(:,1))==1); 
indx2=find(isnan(plot_mat2(:,2))==1); 
indx3=find(isnan(plot_mat2(:,3))==1); 

plot_mat2(indx,1)=0; 
plot_mat2(indx2,2)=0; 
plot_mat2(indx3,3)=0;

[r,c]=size(raw); 
%percent_changed=zeros(5,1); 
group_percent=zeros(4,3); 
group_percent_early=zeros(4,9); 
group_percent_late=zeros(4,9); 
group_total=zeros(4,1); 
group_percent_early_late=zeros(4,9); 
group_stats=zeros(4,3); 

count=1; 
for j=1:r
    spec_num1=g_names{j,2}; 
    spec_1=cell2mat(species(:,1)); 
    indx=find(spec_1==spec_num1); 
    if numel(indx)==1
        spec_num2=species(indx,3); 
        group=cell2mat(spec_num2) ; 
        error=0; 
    else 
        error=1; 
    end         
    if error==0 && group~=0  
        x1=raw(j,2); 
        x2=raw(j,1); 
        diff1=x1/x2;
        if diff1>1.5 || diff1<0.66667
            count=count+1; 
            group_percent(group,1)=group_percent(group,1)+1;
        end 
        x1=raw(j,3); 
        x2=raw(j,1); 
        diff1=x1/x2 ;
        group_total(group)=group_total(group)+1;
        if diff1>1.5 || diff1<0.66667
            count=count+1; 
            group_percent(group,2)=group_percent(group,2)+1;
        end 
        x1=raw(j,2); 
        x2=raw(j,3); 
        diff1=x1/x2; 
        if diff1>1.5 || diff1<0.66667
            count=count+1; 
            group_percent(group,3)=group_percent(group,3)+1;
        end 
        %for each of the 9 samples
        for b=1:9
            x1=EARLY(j,b); 
            x2=LATE(j,b); 
            x3=mean(CONTROL(j,:)); 
            diff1=x1/x3; 
            diff3=x2/x3;
            if diff1>1.5 || diff1<0.6667
                group_percent_early(group,b)=group_percent_early(group,b)+1;
            end 
            if diff3>1.5 || diff3<0.6667
                group_percent_late(group,b)=group_percent_late(group,b)+1; 
            end 
            diff2=x1/x2; 
            if diff2>1.5 || diff2<0.6667
                group_percent_early_late(group,b)=group_percent_early_late(group,b)+1; 
            end 
        end 
    else 
        name; 
    end 
end 
   
bar_mat=zeros(4,3); 
bar_error=zeros(4,3); 
bar_mat2=zeros(4,1); 
bar_error2=zeros(4,1); 
bar_mat3=zeros(3,3); 
bar_error3=zeros(3,3); 
stat_mat=double.empty; 
for j=1:4
    ratio=group_percent_early(j,:)/group_total(j); 
    e_mean=mean(ratio); 
    e_std=std(ratio); 
    bar_mat(j,1)=e_mean; 
    bar_error(j,1)=e_std; 
    group_e=[mean(ratio(1:3)) mean(ratio(4:6)) mean(ratio(7:9))]; 
    ratio=group_percent_late(j,:)/group_total(j) ;
    l_mean=mean(ratio); 
    l_std=std(ratio); 
    bar_mat(j,2)=l_mean; 
    bar_error(j,2)=l_std; 
    group_l=[mean(ratio(1:3)) mean(ratio(4:6)) mean(ratio(7:9))]; 
    ratio=group_percent_early_late(j,:)/group_total(j); 
    b_mean=mean(ratio); 
    b_std=std(ratio); 
    bar_mat(j,3)=b_mean; 
    bar_mat2(j,1)=b_mean; 
    bar_error(j,3)=b_std; 
    bar_error2(j,1)=b_std; 
    group_el=[mean(ratio(1:3)) mean(ratio(4:6)) mean(ratio(7:9))]; 
    [h, stat_mat(j,1)]=ttest2(group_l, group_el); 
    [h, stat_mat(j,2)]=ttest2(group_e, group_el); 
    [h, stat_mat(j,3)]=ttest2(group_e, group_l); 
end 

stat_mat2=double.empty; 
for j=1:3
    if (j<3)
        ratio=group_percent_early(j,:)/group_total(j); 
        e_mean=mean(ratio); 
        e_std=std(ratio); 
        bar_mat3(j,1)=e_mean; 
        bar_error3(j,1)=e_std; 
        group_e=[mean(ratio(1:3)) mean(ratio(4:6)) mean(ratio(7:9))]; 
        ratio=group_percent_late(j,:)/group_total(j) ;
        l_mean=mean(ratio); 
        l_std=std(ratio); 
        bar_mat3(j,2)=l_mean; 
        bar_error3(j,2)=l_std; 
        group_l=[mean(ratio(1:3)) mean(ratio(4:6)) mean(ratio(7:9))]; 
        ratio=group_percent_early_late(j,:)/group_total(j); 
        b_mean=mean(ratio); 
        b_std=std(ratio); 
        bar_mat3(j,3)=b_mean; 
        bar_error3(j,3)=b_std; 
        group_el=[mean(ratio(1:3)) mean(ratio(4:6)) mean(ratio(7:9))]; 
        [h, stat_mat2(j,1)]=ttest2(group_l, group_el); 
        [h, stat_mat2(j,2)]=ttest2(group_e, group_el); 
        [h, stat_mat2(j,3)]=ttest2(group_e, group_l); 
    else
        ge=group_percent_early(3,:)+group_percent_early(4,:);
        gt=group_total(3)+group_total(4); 
        gl=group_percent_late(3,:)+group_percent_late(4,:); 
        gel=group_percent_early_late(3,:)+group_percent_early_late(4,:); 
        
        ratio=ge/gt; 
        e_mean=mean(ratio); 
        e_std=std(ratio); 
        bar_mat3(j,1)=e_mean; 
        bar_error3(j,1)=e_std; 
        group_e=[mean(ratio(1:3)) mean(ratio(4:6)) mean(ratio(7:9))]; 
        ratio=gl/gt ;
        l_mean=mean(ratio); 
        l_std=std(ratio); 
        bar_mat3(j,2)=l_mean; 
        bar_error3(j,2)=l_std; 
        group_l=[mean(ratio(1:3)) mean(ratio(4:6)) mean(ratio(7:9))]; 
        ratio=gel/gt; 
        b_mean=mean(ratio); 
        b_std=std(ratio); 
        bar_mat3(j,3)=b_mean; 
        bar_error3(j,3)=b_std; 
        group_el=[mean(ratio(1:3)) mean(ratio(4:6)) mean(ratio(7:9))]; 
        [h, stat_mat2(j,1)]=ttest2(group_l, group_el); 
        [h, stat_mat2(j,2)]=ttest2(group_e, group_el); 
        [h, stat_mat2(j,3)]=ttest2(group_e, group_l); 
        
    end 
end 

%%%TG and DG plots--------------------------------------------------
%just early-late comparison
bar_mat_=bar_mat2*-100; 
error_mat=bar_error2*100; 
bh=bar(bar_mat_, 1); 
bar_color=[-1 0.25 1 -0.35];
ch=get(bh,'Children');  
set (ch,'CData',bar_color) ; 
colormap (gray); 
hold on 
errorbar (bar_mat_, error_mat, '.k', 'MarkerSize', 2);
hold off 
colormap (gray); 
xlim([0 5]); 
%legend('PL', 'SL', 'NL', 'Location', 'SouthWest'); 
set (gca, 'XTickLabel', {'PL', 'SL', 'TG', 'DG'}); 
set (gca, 'YTick', [-100 -90 -80 -70 -60 -50 -40 -30 -20 -10 0]); 
set (gca, 'YTickLabel', [100 90 80 70 60 50 40 30 20 10 0]); 
ylabel('% of lipid class for species with >1.5X change'); 
set(gca,'FontSize',16,'linewidth',2)
print (gcf, '-dpng','gametocyte_heatmap_bars.png'); 
print (gcf, '-depsc2', ['gametocyte_heatmap_bars.eps']); 
%saveas (gcf, 'gametocyte_heatmap_bars.fig'); 
close all; 

%early-late-control comparison
bar_mat_=bar_mat*-100; 
error_mat=bar_error*100; 
x=1:4; 
for j=1:3
    bar_mat_(:,j); 
    bh=bar(x, bar_mat_(:,j), 1); 
    ch=get(bh,'Children');  
    set (ch,'CData',bar_color) ; 
    colormap (gray); 
    hold on 
    errorbar (x, bar_mat_(:,j), error_mat (:,j), '.k', 'MarkerSize', 2);
	colormap (gray); 
    x=x+[5 5 5 5] ;
end  
hold off
set(gca, 'XTick', [2.5, 7.5, 12.5]); 
set (gca, 'XTickLabel', {'early-control', 'late-control', 'early-late'}); 
set (gca, 'YTick', [-100 -90 -80 -70 -60 -50 -40 -30 -20 -10 0]); 
set (gca, 'YTickLabel', [100 90 80 70 60 50 40 30 20 10 0]); 
set(gca,'FontSize',16,'linewidth',2)
ylabel('% of lipid class for species with >1.5X change'); 
print (gcf, '-dpng','gametocyte_heatmap_control_bars.png'); 
print (gcf, '-depsc2', ['gametocyte_heatmap_control_bars.eps']); 
%saveas (gcf, 'gametocyte_heatmap_control_bars.fig'); 

%NL plots----------------------------------------------

%just early-late comparison
%bar_mat_=bar_mat3(:,3)*-100 ;
bar_mat_=bar_mat2(1:3,:)*-100; 
error_mat=bar_error2(1:3,:)*100; 
bh=bar(bar_mat_, 1); 
bar_color=[-1 0.25 1];
ch=get(bh,'Children');  
set (ch,'CData',bar_color) ; 
colormap (gray); 
hold on 
errorbar (bar_mat_, error_mat, '.k', 'MarkerSize', 2);
hold off 
colormap (gray); 
xlim([0 4]); 
set (gca, 'XTickLabel', {'PL', 'SL', 'NL'}); 
set (gca, 'YTick', [-100 -90 -80 -70 -60 -50 -40 -30 -20 -10 0]); 
set (gca, 'YTickLabel', [100 90 80 70 60 50 40 30 20 10 0]); 
set(gca,'FontSize',16,'linewidth',2)
ylabel('% of lipid class for species with >1.5X change'); 
print (gcf, '-dpng','gametocyte_heatmap_bars_NL.png'); 
print (gcf, '-depsc2', ['gametocyte_heatmap_bars_NL.eps']); 
%saveas (gcf, 'gametocyte_heatmap_bars_NL.fig'); 
close all; 

%early-late-control comparison
%bar_mat_=bar_mat3*-100; 
%error_mat=bar_error3*100; 
bar_mat_=bar_mat(1:3,:)*-100; 
error_mat=bar_error(1:3,:)*100; 
x=1:3; 
for j=1:3
    bh=bar(x, bar_mat_(:,j), 1); 
    ch=get(bh,'Children');  
    set (ch,'CData',bar_color) ; 
    colormap (gray); 
    hold on 
    errorbar (x, bar_mat_(:,j), error_mat (:,j), '.k', 'MarkerSize', 2);
	colormap (gray); 
    x=x+[4 4 4] ;
end  
hold off
set(gca, 'XTick', [1.5, 5.5, 9.5]); 
set (gca, 'XTickLabel', {'early-control', 'late-control', 'early-late'}); 
set (gca, 'YTick', [-100 -90 -80 -70 -60 -50 -40 -30 -20 -10 0]); 
set (gca, 'YTickLabel', [100 90 80 70 60 50 40 30 20 10 0]); 
ylabel('% of lipid class for species with >1.5X change'); 
set(gca,'FontSize',16,'linewidth',2)
print (gcf, '-dpng','gametocyte_heatmap_control_bars_NL.png'); 
print (gcf, '-depsc2',  'gametocyte_heatmap_control_bars_NL.eps'); 
%saveas (gcf, 'gametocyte_heatmap_control_bars_NL.fig'); 
A=1; 

end
