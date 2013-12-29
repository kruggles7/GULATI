function [ plot_mat, bar_mat, percent_changed, G_ ] = lipid_heatmap( mat,names,label, species )

%from the drosophila paper the lipid levels are 
% Mol%= the lipid1(mol)/sum of all membrane lipid species (mol)
%Mol% of lipid class = sum of lipid spces in class/sum of all lipid species
%Mol%moving average = average mol%of lipid species across time
%deviation from moving average=log10(mol%time/mol%movingaverage)
%  mat=super_;
%  names=super_names;
%  label='super'; 

G=cell.empty; 
group={'PL', 'SL', 'NL'}; 
match_=cell2mat(species(:,3));
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
        if numel(indx)==1
            group_=species{indx,3};   %%group they belong to  
            %group_=str2num(group_); 
            x1=plot_mat(j,k); 
            x2=plot_mat(j,k+1); 
            diff1=x1/x2; 
            diff2=x1/x2; 
            if diff1>1.5 || diff2>1.5
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
                diff1=x1/x2; 
                diff2=x1/x2;
                if diff1>1.5 || diff2>1.5 
                    group_percent_all(group_,k,b)=group_percent_all(group_,k,b)+1;
                end 
            end 
                
            %group1
            x1=plot_mat1(j,k); 
            x2=plot_mat1(j,k+1); 
            diff1=x1/x2; 
            diff2=x1/x2;
            if diff1>1.5 || diff2>1.5
                count=count+1; 
                group_percent1(group_,k)=group_percent1(group_,k)+1;
            end 
            
            %group2
            x1=plot_mat2(j,k); 
            x2=plot_mat2(j,k+1); 
            diff1=x1/x2; 
            diff2=x1/x2;
            if diff1>1.5 || diff2>1.5
                count=count+1; 
                group_percent2(group_,k)=group_percent2(group_,k)+1;
            end
            
            %group3
            x1=plot_mat3(j,k); 
            x2=plot_mat3(j,k+1); 
            diff1=x1/x2; 
            diff2=x1/x2;
            if diff1>1.5 || diff2>1.5
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
        bar_mat1(c_b)=group_percent1(j,k)/group_total(j)*100; 
        bar_mat2(c_b)=group_percent2(j,k)/group_total(j)*100; 
        bar_mat3(c_b)=group_percent3(j,k)/group_total(j)*100; 
        for b=1:9
            bar_mat_all(c_b,b)=group_percent_all(j,k,b)/group_total(j)*100;
        end 
        test=[bar_mat1(c_b) bar_mat2(c_b) bar_mat3(c_b)]; 
        bar_avg(c_b)=mean(test); 
        bar_std(c_b)=std(test); 
        test2=[bar_mat_all(c_b,1) bar_mat_all(c_b,2) bar_mat_all(c_b,3) bar_mat_all(c_b,4) bar_mat_all(c_b,5) bar_mat_all(c_b,6) bar_mat_all(c_b,7) bar_mat_all(c_b,8) bar_mat_all(c_b,9)]; 
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
        count=x; 
        for j=1:5
            group2=bar_mat_all(count,:); 
            [h,p]=ttest(group1, group2); 
            bar_all_sig(k,j)=p; 
            count=count+4;
        end 
        x=x+1; 
    end 
    count2=count2+4; 
end 
     

%bar_mat_=rot90(bar_mat); 
%bar_mat_=fliplr(bar_mat_) ;
bar_mat=bar_mat*-1; 
bh=bar(bar_mat, 1); 
bar_color=[-1 0 1 0 -1 0 1 0 -1 0 1 0 -1 0 1 0 -1 0 1 0];
ch=get(bh,'Children');  
set (ch,'CData',bar_color) ; 
colormap (gray); 
xlim([0 21]); 
%legend('PL', 'SL', 'NL', 'Location', 'SouthWest'); 
set (gca, 'XTickLabel', []); 
set (gca, 'YTick', [-50 -40 -30 -20 -10 0]); 
set (gca, 'YTickLabel', [50 40 30 20 10 0]); 
ylabel('% of lipid class for species with >1.5X change'); 
print (gcf, '-dpng',[label '_heatmap_bars.png']); 

%AVG of triplicates
bar_avg=bar_avg*-1; 
bh=bar(bar_avg, 1); 
bar_color=[-1 0 1 0 -1 0 1 0 -1 0 1 0 -1 0 1 0 -1 0 1 0];
ch=get(bh,'Children');  
set (ch,'CData',bar_color) ; 
colormap (gray); 
hold on 
errorbar (bar_avg, bar_std, '.k', 'MarkerSize', 2);
hold off
xlim([0 21]); 
%legend('PL', 'SL', 'NL', 'Location', 'SouthWest'); 
set (gca, 'XTickLabel', []); 
set (gca, 'YTick', [-50 -40 -30 -20 -10 0]); 
set (gca, 'YTickLabel', [50 40 30 20 10 0]); 
ylabel('% of lipid class for species with >1.5X change'); 
print (gcf, '-dpng',[label '_heatmap_bars_triplicate.png']); 

%AVG of 9
bar_avg_all=bar_avg_all*-1; 
bh=bar(bar_avg_all, 1); 
bar_color=[-1 0 1 0 -1 0 1 0 -1 0 1 0 -1 0 1 0 -1 0 1 0];
ch=get(bh,'Children');  
set (ch,'CData',bar_color) ; 
colormap (gray); 
hold on 
errorbar (bar_avg_all, bar_std_all, '.k', 'MarkerSize', 2);
hold off
xlim([0 21]); 
%legend('PL', 'SL', 'NL', 'Location', 'SouthWest'); 
set (gca, 'XTickLabel', []); 
set (gca, 'YTick', [-50 -40 -30 -20 -10 0]); 
set (gca, 'YTickLabel', [50 40 30 20 10 0]); 
ylabel('% of lipid class for species with >1.5X change'); 
print (gcf, '-dpng',[label '_heatmap_bars_all.png']); 

%AVG of 9 no bars
bh=bar(bar_avg_all, 1); 
bar_color=[-1 0 1 0 -1 0 1 0 -1 0 1 0 -1 0 1 0 -1 0 1 0];
ch=get(bh,'Children');  
set (ch,'CData',bar_color) ; 
colormap (gray); 
xlim([0 21]); 
%legend('PL', 'SL', 'NL', 'Location', 'SouthWest'); 
set (gca, 'XTickLabel', []); 
set (gca, 'YTick', [-50 -40 -30 -20 -10 0]); 
set (gca, 'YTickLabel', [50 40 30 20 10 0]); 
ylabel('% of lipid class for species with >1.5X change'); 
print (gcf, '-dpng',[label '_heatmap_bars_all_nobars.png']); 


indx=find(plot_mat==0); 
plot_mat(indx)=1; 
indx=find(isnan(plot_mat)==1); 
plot_mat(indx)=1; 
plot_mat=log10(plot_mat); 
x={'8','16','24','32','40','48'}; 
plots=flipud(plot_mat); 
G=flipud(G); 
name_group=flipud(names(:,1)); 
cmap=redgreencmap(256); 
H=HeatMap(plots,'ColumnLabels',x,'RowLabels',name_group, 'DisplayRange', 3, 'Colormap', cmap, 'Symmetric','true'); 
%H=heatmap_rg(plots, x, name_group, [],-1, 1, 'ColorBar', 1, 'Colormap', 'money');%'ColorLevels', 256, 
%H=heatmap(plots, x, name_group, [], 'ColorBar', 1, 'Colormap', 'money');%'ColorLevels', 256, 
plot(H);
print (gcf, '-dpng',[ label '_heatmap.png']); 

%sorted-------------------------------------------------------------------
[r,c]=size(plots); 
S(:,1)=1:r; 
S(:,2)=plots(:,1)-plots(:,6); 
S=sortrows(S,2); 
plots2=zeros(r,c); 
names2=cell.empty; 
G_=double.empty; 
for i=1:r
    indx=S(i,1); 
    plots2(i,:)=plots(indx,:);
    temp=G{indx,1}; 
    G_(i,1)=temp; 
end 
H=HeatMap(plots2,'ColumnLabels',x, 'RowLabels', G_, 'DisplayRange', 3, 'Colormap', cmap, 'Symmetric','true');  
plot(H);
print (gcf, '-dpng',[ label '_sorted_heatmap.png']); 

%clustergram(plots, 'Cluster', 1); 
close all 
%G

% max_=max(cell2mat(names(:,2))); 
% count=1; 
% for i=1:max_
%     indx=find(cell2mat(names(:,2))==i); 
%     plot_mat_=plot_mat(indx,:); 
%     plot_mat_=flipud(plot_mat_); 
%     name_=names(indx,1); 
%     name_=flipud(name_); 
%     H=HeatMap(plot_mat_,'ColumnLabels',x,'RowLabels',name_, 'DisplayRange',1); 
%     plot(H); 
%     i_str=num2str(i); 
%     print (gcf, '-dpng',[i_str '_' label '_heatmap_species.png']); 
%     count=count+1; 
%     %close all 
% end 
end

