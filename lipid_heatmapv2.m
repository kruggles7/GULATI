function [ plot_mat, bar_mat, stat_final, number_changed, bar_all_sig, names_] = lipid_heatmapv2( mat,names,label,species)

%from the drosophila paper the lipid levels are 
% Mol%= the lipid1(mol)/sum of all membrane lipid species (mol)
%Mol% of lipid class = sum of lipid spces in class/sum of all lipid species
%Mol%moving average = average mol%of lipid species across time
%deviation from moving average=log10(mol%time/mol%movingaverage)
%   mat=DD2f_R_;
%   names=DD2f_names;
%   label='RBC'; 

stat_final=double.empty; 
GRAY=[190/256 190/256 190/256]; 
G=cell.empty; 
group=species(:,2); 
match_=cell2mat(species(:,1)); 
[r,c]=size(names); 
names_=zeros(r,1); 

[r,c]=size(species); 
group=cell2mat(species(:,1)) 
group_=cell2mat(species(:,3)); 
n=cell2mat(names(:,2)); 
for i=1:r
    species_group=group(i); 
    indx=find(n==species_group); 
    names_(indx)=group_(i); 
end 

number_changed=double.empty; 




group_percent=zeros(4,5); 
group_total=zeros(4,1); 
group_percent1=zeros(4,5); 
group_percent2=zeros(4,5); 
group_percent3=zeros(4,5); 
group_percent_all=zeros(4,5,9); 
[r,c]=size(mat); 
n=c/6; 
mol_time=zeros(r,c); 
sum_mat=sum(mat); 
%name_group=cell2mat(names(:,2)); 
%unique_group=unique(name_group); 
%find the mol%time for each lipid 
for k=1:r
    for j=1:c
        mol_time(k,j)=mat(k,j)/sum_mat(1,j); 
    end 
end

mov_avg=mean(mol_time,2); 

%deviation from moving average
for k=1:r
    for j=1:c
        dev_moving(k,j)=mol_time(k,j)/mov_avg(k); 
    end 
end

%foreach time
plot_mat=double.empty; 

%for each experiment
plot_mat1=double.empty; 
plot_mat2=double.empty; 
plot_mat3=double.empty; 
i=1; 
n=1; 
plot_mat_all=double.empty; 
for k=1:6 
    for j=1:r
        plot_mat(j,k)=mean(dev_moving(j,i:i+8)); 
    end
    for j=1:r
        plot_mat1(j,k)=mean(dev_moving(j,i:i+2)); 
        plot_mat2(j,k)=mean(dev_moving(j,i+3:i+5));
        plot_mat3(j,k)=mean(dev_moving(j,i+6:i+8)); 
        for n=1:9
            plot_mat_all(j,k,n)=dev_moving(j,i+n-1); 
        end 
    end 
    i=i+9; 
end

% %statistics on the dev_moving 
% stat_final=cell(r,3); 
% for j=1:r %each row
%     stat_mat=double.empty; 
%     stat_final{j,1}=names{j,1}; 
%     i=1; 
%     for k=1:6
%         stat_mat(:,k)=dev_moving(j,i:i+8)'; 
%         i=i+9; 
%     end 
%     indx=find(isnan(stat_mat)==0); 
%     if numel(indx>0)
%         [p,a,s]=anova1(stat_mat); 
%         close all
%         stat_final{j,2}=p; 
%         %ttest
%         temp_ttest=[]; 
%         for k=1:6
%             group1=stat_mat(:,k); 
%             for i=(k+1):6
%                 group2=stat_mat(:,i); 
%                 [hval, pval]=ttest2(group1, group2); 
%                 if pval<0.05
%                     temp1=num2str(k); 
%                     temp2=num2str(i); 
%                     temp_ttest=[temp_ttest temp1 '/' temp2 ', ']; 
%                 end 
%             end 
%         end 
%         stat_final{j,3}=temp_ttest; 
%     else 
%         stat_final{j,2}=NaN; 
%         stat_final{j,3}=NaN; 
%     end 
% end 

%find the log ratio ! 
indx=find(plot_mat==0); 
plot_mat(indx)=1; 
indx=find(isnan(plot_mat)==1); 
plot_mat(indx)=1; 
plot_mat=log2(plot_mat); 

indx=find(plot_mat_all==0); 
plot_mat_all(indx)=1; 
indx=find(isnan(plot_mat_all)==1); 
plot_mat_all(indx)=1; 
plot_mat_all=log2(plot_mat_all); 

indx=find(plot_mat1==0); 
plot_mat1(indx)=1; 
indx=find(isnan(plot_mat1)==1); 
plot_mat1(indx)=1; 
plot_mat1=log2(plot_mat1); 

indx=find(plot_mat2==0); 
plot_mat2(indx)=1; 
indx=find(isnan(plot_mat2)==1); 
plot_mat2(indx)=1; 
plot_mat2=log2(plot_mat2); 

indx=find(plot_mat3==0); 
plot_mat3(indx)=1; 
indx=find(isnan(plot_mat3)==1); 
plot_mat3(indx)=1; 
plot_mat3=log2(plot_mat3);

%check for number of PL that have > 1.5x change 
percent_changed=zeros(5,1);

for k=1:5 
    count=1; 
    for j=1:r
        %FIND THE NAME OF EACH 
        text_=names{j,1}; 
        name=''; 
        k_=strfind(text_,' '); 
        k1=(k_(1)-1); 
        for p_=1:k1
            c_=text_(p_);
            name=[name c_];
        end 
        if strcmp (name, 'Acyl')==1 || strcmp(name,'Lyso')==1
            for m_=p_+1:p_+3
                c_=text_(m_); 
                name=[name c_]; 
            end 
        end 
        indx=find(strcmp(name,species(:,2))==1); 
        G{j,1}=species{indx,3};
        H2{j,1}=species{indx,1}; 
        if numel(indx)==1
            group_=species{indx,3};   %%group they belong to  
            group2=species{indx,1}; 
            %group_=str2num(group_); 
            x1=plot_mat(j,k); 
            x2=plot_mat(j,k+1); 
            diff1=abs(x2-x1); 
            if diff1>log2(1.5)
                count=count+1; 
                group_percent(group_,k)=group_percent(group_,k)+1;
            end 
            if k==1
                group_total(group_)=group_total(group_)+1; 
            end 
            
            %all  plot_mat_all(j,k,i)
            for b=1:9
                x1=plot_mat_all(j,k,b); 
                x2=plot_mat_all(j,k+1,b); 
                diff1=abs(x2-x1); 
                if diff1>log2(1.5)
                    group_percent_all(group_,k,b)=group_percent_all(group_,k,b)+1; %group_percent_all(GROUP, TIME, N)
                end 
            end 
                
            %group1
            x1=plot_mat1(j,k); 
            x2=plot_mat1(j,k+1); 
            diff1=abs(x2-x1); 
            if diff1>log2(1.5)
                count=count+1; 
                group_percent1(group_,k)=group_percent1(group_,k)+1;
            end 
            
            %group2
            x1=plot_mat2(j,k); 
            x2=plot_mat2(j,k+1); 
            diff1=abs(x2-x1); 
            if diff1>log2(1.5)
                count=count+1; 
                group_percent2(group_,k)=group_percent2(group_,k)+1;
            end
            
            %group3
            x1=plot_mat3(j,k); 
            x2=plot_mat3(j,k+1); 
            diff1=abs(x2-x1); 
            if diff1>log2(1.5)
                count=count+1; 
                group_percent3(group_,k)=group_percent3(group_,k)+1;
            end  
        else 
            name; 
        end 
    end 
end 
           
            
%MAKE the >1.5 plot
bar_mat=zeros(20,1); 
bar_mat1=zeros(20,1); 
bar_mat2=zeros(20,1); 
bar_mat3=zeros(20,1); 
bar_mat_all=zeros(20,9); 
bar_sig=zeros(20,1); 
bar_avg=zeros(20,1); 
bar_std=zeros(20,1); 
bar_avg_all=zeros(20,1); 
bar_std_all=zeros(20,1); 
c_b=1; 
for k=1:5
    for j=1:3
        bar_mat(c_b)=group_percent(j,k)/group_total(j)*100; 
        number_changed(c_b)=group_percent(j,k); 
        bar_mat1(c_b)=group_percent1(j,k)/group_total(j)*100; 
        bar_mat2(c_b)=group_percent2(j,k)/group_total(j)*100; 
        bar_mat3(c_b)=group_percent3(j,k)/group_total(j)*100; 
        for b=1:9
            bar_mat_all(c_b,b)=group_percent_all(j,k,b)/group_total(j)*100; %group_percent_all(GROUP, TIME, N)
            %bar_mat_all(GROUP, SAMPLE) 1-5, skip, 7-11, skip, 13-17
        end 
        test=[bar_mat1(c_b) bar_mat2(c_b) bar_mat3(c_b)]; 
        bar_avg(c_b)=mean(test); 
        bar_std(c_b)=std(test); 
        test2=bar_mat_all(c_b,:); 
        bar_avg_all(c_b)=mean(test2); 
        bar_std_all(c_b)=std(test2); 
        c_b=c_b+1; 
    end 
    bar_mat(c_b)=0; 
    bar_mat1(c_b)=0; 
    bar_mat2(c_b)=0; 
    bar_mat3(c_b)=0; 
    bar_std(c_b)=0; 
    bar_avg(c_b)=0; 
    bar_std_all(c_b)=0; 
    bar_avg_all(c_b)=0; 
    c_b=c_b+1; 
end 

%test significance across experiments
bar_all_sig=double.empty; 
count2=1; 
for i=1:5
    x=1; 
    for k=count2:count2+2
        group1=bar_mat_all(k,:);  
        count=x; %x= 1,2,3,4,5 
        for j=1:5
            group2=bar_mat_all(count,:);  %bar_mat_all(GROUP, SAMPLE) 1-5, skip, 7-11, skip, 13-17
            [h,p]=ttest2(group1, group2); 
            bar_all_sig(k,j)=p; %for each of the 5 groups 
            count=count+4; %if count=1
        end 
        x=x+1; 
    end 
    count2=count2+4; %compare every 5 to make sure you compare the correct combiations
end 


cd ..
cd results
cd heatmaps

%calculated together----------------------------------------------
%bar_mat_=rot90(bar_mat); 
%bar_mat_=fliplr(bar_mat_) ;
bar_mat=bar_mat*-1; 
x=1:4:17; 
for i=1:numel(x)
    bar(x(i),bar_mat(x(i)),1,'FaceColor','k', 'EdgeColor', 'k'); 
	hold on
end
x=2:4:18; 
for i=1:numel(x)
    bar(x(i),bar_mat(x(i)),1,'FaceColor', GRAY, 'EdgeColor', 'k'); 
end 
x=3:4:19; 
for i=1:numel(x)
    bar(x(i),bar_mat(x(i)),1,'FaceColor', 'w', 'EdgeColor', 'k'); 
end
hold off
xlim([0 21]); 
%legend('PL', 'SL', 'NL', 'Location', 'SouthWest'); 
set (gca, 'XTickLabel', []); 
set (gca, 'YTick',[-100 -80 -60 -40 -20 0]); 
set (gca, 'YTickLabel', [100 80 60 40 20 0]); 
ylim([-120 0])
ylabel('% of lipid class for species with >1.5X change'); 
set(gca,'FontSize',16,'linewidth',2)
print (gcf, '-dpng',[label '_heatmap_bars.png']); 
print (gcf, '-depsc2', [label '_heatmap_bars.eps']); 
%saveas (gcf,[ label '_heatmap_bars.fig']); 

%AVG of triplicates----------------------------------------------------
bar_avg=bar_avg*-1; 
x=1:4:17; 
for i=1:numel(x)
    bar(x(i),bar_avg(x(i)),1,'FaceColor','k', 'EdgeColor', 'k'); 
	hold on
    errorbar (x(i), bar_avg(x(i)), bar_std(x(i)), '.k', 'MarkerSize', 2,'linewidth',1);
end
x=2:4:18; 
for i=1:numel(x)
    bar(x(i),bar_avg(x(i)),1,'FaceColor', GRAY, 'EdgeColor', 'k'); 
    errorbar (x(i), bar_avg(x(i)), bar_std(x(i)), '.k', 'MarkerSize', 2,'linewidth',1);
end 
x=3:4:19; 
for i=1:numel(x)
    bar(x(i),bar_avg(x(i)),1,'FaceColor', 'w', 'EdgeColor', 'k'); 
    errorbar (x(i), bar_avg(x(i)), bar_std(x(i)), '.k', 'MarkerSize', 2,'linewidth',1);
end
hold off
xlim([0 21]); 
%legend('PL', 'SL', 'NL', 'Location', 'SouthWest'); 
set (gca, 'XTickLabel', []); 
set (gca, 'YTick',[-100 -80 -60 -40 -20 0]); 
set (gca, 'YTickLabel', [100 80 60 40 20 0]); 
ylim ([-120 0]); 
ylabel('% of lipid class for species with >1.5X change'); 
set(gca,'FontSize',16,'linewidth',2)
print (gcf, '-dpng',[label '_heatmap_bars_triplicate.png']); 
print (gcf, '-depsc2', [label '_heatmap_bars_triplicate.eps']); 
%saveas (gcf,[ label '_heatmap_bars_triplicate.fig']); 

%AVG of 9 with bars--------------------------------------------------------
bar_avg_all=bar_avg_all*-1; 
%AVG of 9 no bars
x=1:4:17; 
for i=1:numel(x)
    bar(x(i),bar_avg_all(x(i)),1,'FaceColor','k', 'EdgeColor', 'k'); 
	hold on
    errorbar (x(i), bar_avg_all(x(i)), bar_std_all(x(i)), '.k', 'MarkerSize', 2,'linewidth',1);
end
x=2:4:18; 
for i=1:numel(x)
    bar(x(i),bar_avg_all(x(i)),1,'FaceColor', GRAY, 'EdgeColor', 'k'); 
    errorbar (x(i), bar_avg_all(x(i)), bar_std_all(x(i)), '.k', 'MarkerSize', 2,'linewidth',1);
end 
x=3:4:19; 
for i=1:numel(x)
    bar(x(i),bar_avg_all(x(i)),1,'FaceColor', 'w', 'EdgeColor', 'k'); 
    errorbar (x(i), bar_avg_all(x(i)), bar_std_all(x(i)), '.k', 'MarkerSize', 2,'linewidth',1);
end
hold off
xlim([0 21]); 

%legend('PL', 'SL', 'NL', 'Location', 'SouthWest'); 
set (gca, 'XTickLabel', []); 
set (gca, 'YTick',[-100 -80 -60 -40 -20 0]); 
set (gca, 'YTickLabel', [100 80 60 40 20 0]); 
ylim([-120 0])
ylabel('% of lipid class for species with >1.5X change'); 
set(gca,'FontSize',16,'linewidth',2)
print (gcf, '-dpng',[label '_heatmap_bars_all.png']); 
print (gcf, '-depsc2', [label '_heatmap_bars_all.eps']); 
%saveas (gcf,[ label '_heatmap_bars_all.fig']); 

%AVG of 9 no bars---------------------------------------------------------
x=1:4:17; 
for i=1:numel(x)
    bar(x(i),bar_avg_all(x(i)),1,'FaceColor','k', 'EdgeColor', 'k'); 
	hold on
end
x=2:4:18; 
for i=1:numel(x)
    bar(x(i),bar_avg_all(x(i)),1,'FaceColor', GRAY, 'EdgeColor', 'k'); 
end 
x=3:4:19; 
for i=1:numel(x)
    bar(x(i),bar_avg_all(x(i)),1,'FaceColor', 'w', 'EdgeColor', 'k'); 
end
hold off
colormap (gray); 
xlim([0 21]); 
%legend('PL', 'SL', 'NL', 'Location', 'SouthWest'); 
set (gca, 'XTickLabel', []); 
set (gca, 'YTick',[-100 -80 -60 -40 -20 0]); 
set (gca, 'YTickLabel', [100 80 60 40 20 0]); 
ylabel('% of lipid class for species with >1.5X change'); 
set(gca,'FontSize',16,'linewidth',2)
ylim([-120 0])
print (gcf, '-dpng',[label '_heatmap_bars_all_nobars.png']); 
print (gcf, '-depsc2', [label '_heatmap_bars_all_nobars.eps']); 
%saveas (gcf,[ label '_heatmap_bars_all_nobars.fig']); 

%HEATMAP---------------------------------------------------------------

x={'8','16','24','32','40','48'}; 
plots=flipud(plot_mat); 
%G=flipud(G); 
name_group=names(:,1); 
cmap=redgreencmap(256); 
%cmap='jet'; 
H=HeatMap(plots,'ColumnLabels',x,'RowLabels',name_group, 'DisplayRange', 4, 'Colormap', cmap, 'Symmetric','true'); 
set(gca,'FontSize',16,'linewidth',2)
plot(H);
ylim([-100 0])
print (gcf, '-dpng',[ label '_heatmap.png']); 
print (gcf, '-depsc2', [label '_heatmap.eps']); 
saveas(gcf, [label '_heatmap.fig']); 


cd .. 
cd ..
cd programs

end

