function [anova_sat, anova_chain, anova_combo, plot_combo, combo_text] = sat_length_gamsv3( mat, lipids, species, label, T )
%UNTITLED Summary of this function goes here
% %   Detailed explanation goes here
% 
% mat=all; 
% lipids=gam_text; 
% label='gams'; 
% T=3; 
[r,c]=size(mat);  

[r,c]=size(lipids); 
group={'PL', 'SL', 'NL'}; 
match=cell2mat(species(:,3)); 
anova_sat=cell.empty;
anova_chain=cell.empty; 
anova_combo=cell.empty; 
a_count=1; 
c_count=1; 
co_count=1; 
cd ..
cd results
mkdir saturation
cd saturation

match_=match; 
match_(match==4)=3; 
for i=1:3
    indx=find(match_==i); 
    L=species(indx,2); 
    saturation=double.empty; 
    chain_length=double.empty;
    time_mat=double.empty; 
    col_label=double.empty; 
    levels=double.empty; 
    n_=1; 
    for k=1:length(L)
        lipid=L{k}; %the name of the lipid we are including
        [r,c]=size(lipids); 
        for j=1:r
            text_=lipids{j,1}; 
            text1=text_; 
            name=''; 
            k_=strfind(text_,' '); 
            k1=(k_(1)-1); 
            for p_=1:k1
                c_=text_(p_);
                name=[name c_];
            end 
            if strcmp (name, 'Acyl')==1 || strcmp(name, 'Lyso')==1
                for m_=p_+1:p_+3
                    c_=text_(m_); 
                    name=[name c_]; 
                end 
            end 
            q=strfind(text_,':'); 
            if numel(q)>0 && numel(q)<3
                q1=q(1)-1;
            else 
                q1=0; 
            end 
            if strcmp(name,lipid)==1 %same family
                %look at the chain length   
                if q1>0
                    chain=''; 
                    sat_=''; 
                    Q=strfind(text_, '/'); 
                    P=strfind(text_,'-');
                    if numel(Q)>0
                        %q1=q(2)-1; 
                        rL=Q-1; 
                    elseif numel(P)>0
                        %q1=q(2)-1;
                        rL=P(1)-1; 
                    else 
                        rL=length(text_);
                    end 
                    for p_=q1-1:q1
                        c_=text_(p_);  
                        chain=[chain c_] ;
                    end
                    chain = str2num(chain); 
                    for p_=q1+2:rL
                        c_=text_(p_) ;
                        sat_=[sat_ c_]; 
                    end 
                    sat_=str2num(sat_);
                    
               %FOR ADDING SECOND NUMBER -----------------
                    rL=length(text_); 
                    if numel(Q)>0
                        q2=q(2)-1; 
                        %rL=Q-1; 
                    elseif numel(P)>0
                        q2=q(2)-1;
                        %rL=P-1; 
                    else 
                        q2=0; 
                    end 
                    if i==2 %PL
                        chain2=''; 
                        sat2=''; 
                        for p_=q2-1:q2
                            c_=text_(p_);  
                            chain2=[chain2 c_] ;
                        end
                        chain2 = str2num(chain2); 
                        chain_length(n_,1)=chain + chain2; 

                        for p_=q2+2:rL
                            c_=text_(p_) ;
                            sat2=[sat2 c_]; 
                        end 
                        sat2=str2num(sat2);
                        saturation(n_,1)=sat_+sat2;
                    elseif numel(P)>0
                        chain2=''; 
                        sat2=''; 
                        for p_=q2-1:q2
                            c_=text_(p_);  
                            chain2=[chain2 c_] ;
                        end
                        chain2 = str2num(chain2); 
                        chain_length(n_,1)=chain2; 

                        for p_=q2+2:rL
                            c_=text_(p_) ;
                            sat2=[sat2 c_]; 
                        end 
                        sat2=str2num(sat2);
                        saturation(n_,1)=sat2;
                        
                    else
                        chain_length(n_,1)=chain; 
                        saturation(n_,1)=sat_; 
                    end 
                    levels(n_,:)=mat(j,:); 
                    col_label{n_,1}=text1;  
                    n_=n_+1;
                end  
            end 
        end
    end 

    %make an average with time for each member of the group
    [r1,c1]=size(levels); 
    levels_sum=nansum(levels); 
    levels_percent=nan(r1,c1); 
    for k2=1:r1
        for j2=1:c1
            levels_percent(k2,j2)=(levels(k2,j2)/levels_sum(j2))*100; 
        end 
    end 
    if T>1 %% ALL EXCEPT FOR CONTROL 
        for k2=1:r1
            d=1; 
            for j=1:T 
                time_mat(k2,j)=nanmean(levels(k2,d:d+8)); 
                d=d+9;  
            end
        end
        sum_=nansum(time_mat);
        unique_sat=unique(saturation); 
        plot_sat=zeros(length(unique_sat),T); 
        text_sat=double.empty; 
        text_length=double.empty; 
        stat_mat=double.empty; 
        for j=1:length(unique_sat)
            u=unique_sat(j); 
            anova_sat{a_count,1}=group{i};
            anova_sat{a_count,2}=u; 
            indx=find(saturation==u); %find the same saturation level
            lipid_new=time_mat(indx,:);
            [rl,cl]=size(lipid_new); 
            if rl>1
                lipid_sum=nansum(lipid_new); 
            else 
                lipid_sum=lipid_new; 
            end 
            for j2=1:T
                plot_sat(j,j2)=(lipid_sum(j2)/sum_(j2))*100;  
            end 
            text_sat(j)=u; 
            stat_sat=levels_percent(indx,:); 
            stat_mat=nan(9,3); 
            [r3,c3]=size(stat_sat); 
            if r3>1
                stat_sat=nansum(stat_sat); 
            end  
            d=1; 
            for j2=1:T 
                temp=stat_sat(d:d+8); 
                temp(temp==0)=nan;  
                stat_mat(:,j2)=temp; 
                d=d+9; 
            end
            p_sat=anova1(stat_mat); 
            close all
            anova_sat{a_count,3}=p_sat; 
            temp_ttest=[]; 
            for n=1:T
                group1=stat_mat(:,n); 
                for m=(n+1):T
                    group2=stat_mat(:,m);  
                    [hval,pval]=ttest2(group1,group2); 
                    if pval<0.01
                        temp_ttest=[temp_ttest n m]; 
                    end 
                end 
            end 
            L_t=length(temp_ttest); 
            temp3=''; 
            for J=1:2:L_t
                temp1=num2str(temp_ttest(J)); 
                temp2=num2str(temp_ttest(J+1)); 
                temp3=[temp3 temp1 '/' temp2 ',']; 
            end 
            anova_sat{a_count,4}=temp3; 
            if isempty(temp3)==1
                anova_sat{a_count,4}='NA'; 
            end 
            a_count=a_count+1;  
        end 
        unique_length=unique(chain_length); 
        plot_length=zeros(length(unique_length),T); 
        for j=1:length(unique_length) 
            u=unique_length(j); 
            anova_chain{c_count,1}=group{i}; 
            anova_chain{c_count,2}=u;
            indx=find(chain_length==u); %find the same saturation level
            lipid_new=time_mat(indx,:);
            [rl,c1]=size(lipid_new); 
            if rl>1
                lipid_sum=nansum(lipid_new); 
            else 
                lipid_sum=lipid_new; 
            end 
            for j2=1:T
                plot_length(j,j2)=(lipid_sum(j2)/sum_(j2))*100;  
            end 
            text_length(j)=u; 
            stat_chain=levels_percent(indx,:); 
            stat_mat=nan(9,3); 
            [r3,c3]=size(stat_chain); 
            if r3>1
                stat_chain=nansum(stat_chain); 
            end  
            d=1; 
            for j2=1:T 
                temp=stat_chain(d:d+8);
                temp(temp==0)=nan; 
                stat_mat(:,j2)=temp; 
                d=d+9; 
            end
            p_chain=anova1(stat_mat); 
            close all
            anova_chain{c_count,3}=p_chain; 
            temp_ttest=[]; 
            for n=1:T
                group1=stat_mat(:,n); 
                for m=(n+1):T
                    group2=stat_mat(:,m); 
                    [hval,pval]=ttest2(group1,group2); 
                    if pval<0.01
                        temp_ttest=[temp_ttest n m]; 
                    end 
                end 
            end 
            L_t=length(temp_ttest); 
            temp3=''; 
            for J=1:2:L_t
                temp1=num2str(temp_ttest(J)); 
                temp2=num2str(temp_ttest(J+1)); 
                temp3=[temp3 temp1 '/' temp2 ',']; 
            end 
            anova_chain{c_count,4}=temp3; 
            if isempty(temp3)==1
                anova_chain{c_count,4}='NA'; 
            end 
            c_count=c_count+1; 
        end 
        plot_combo=double.empty;
        count=1; 
        for j=1:length(unique_sat)
            us=unique_sat(j); 
            indxs=find(saturation==us); 
            for j2=1:length(unique_length)
                error=0; 
                ul=unique_length(j2); 
                indxl=find(chain_length==ul); 
                indx=intersect(indxl,indxs); 
                if numel(indx)==0 
                    error=1; 
                elseif numel(indx)>1
                    lipid_new=nansum(time_mat(indx,:));
                else 
                    lipid_new=time_mat(indx,:); 
                end 
                if error==0 
                    [rl,c1]=size(lipid_new); 
                    if rl>1
                        lipid_sum=nansum(lipid_new); 
                    else 
                        lipid_sum=lipid_new; 
                    end 
                    for j3=1:T
                        plot_combo(count,j3)=(lipid_sum(j3)/sum_(j3))*100; 
                    end 
                    text_combo(count,1)=us;
                    text_combo(count,2)=ul; 
                    count=count+1; 
                    %STATS COMBO------------------------
                    anova_combo{co_count,1}=group{i}; 
                    anova_combo{co_count,2}=ul;
                    anova_combo{co_count,3}=us; 
                    stat_combo=levels_percent(indx,:); 
                    stat_mat=nan(9,3); 
                    [r3,c3]=size(stat_combo); 
                    if r3>1
                        stat_combo=nansum(stat_combo); 
                    end  
                    d=1; 
                    for j3=1:T 
                        temp=stat_combo(d:d+8); 
                        temp(temp==0)=nan; 
                        stat_mat(:,j3)=temp; 
                        d=d+9; 
                    end
                    p_chain=anova1(stat_mat); 
                    close all
                    anova_combo{co_count,4}=p_chain; 
                    temp_ttest=[]; 
                    for n=1:T
                        group1=stat_mat(:,n); 
                        for m=(n+1):T
                            group2=stat_mat(:,m); 
                            [hval,pval]=ttest2(group1,group2); 
                            if pval<0.01
                                temp_ttest=[temp_ttest n m]; 
                            end 
                        end 
                    end 
                    L_t=length(temp_ttest); 
                    temp3=''; 
                    for J=1:2:L_t
                        temp1=num2str(temp_ttest(J)); 
                        temp2=num2str(temp_ttest(J+1)); 
                        temp3=[temp3 temp1 '/' temp2 ',']; 
                    end 
                    anova_combo{co_count,5}=temp3; 
                    if isempty(temp3)==1
                        anova_combo{co_count,5}='NA'; 
                    end 
                    co_count=co_count+1; 
                end 
            end 
        end 
        [r3,c3]=size(plot_combo); 
        plot_combo=rot90(plot_combo); 
        plot_combo=flipud(plot_combo); 
        NAME=group{i}; 
        bar (1:c3, plot_combo, 'stack');
        hold on 
        combo_text=cell.empty; 
        for i3=1:r3
            t=num2str(text_combo(i3,1)); 
            v=num2str(text_combo(i3,2)); 
            x=[v ':' t]; 
            combo_text{i3}=x;  
        end 
        
        legend(combo_text, 'Location', 'EastOutside');    
        ylim([0 100]); 
        ylabel(['% of ' NAME ' with length:saturation']);    
        if T==6
            set(gca, 'XTickLabel', {'8', '16', '24', '32', '40', '48'}); 
        elseif T==3
            set (gca, 'XTickLabel', {'Control', 'Early', 'Late'});  
        end 
        
        print (gcf, '-dpng', [label '_' NAME '_combo.png']); 
        print (gcf, '-depsc2', [label '_' NAME '_combo.eps']);
        close 

        [r3,c3]=size(plot_length); 
        plot_length=rot90(plot_length); 
        plot_length=flipud(plot_length); 
        bar (1:c3, plot_length, 'stack'); 
        hold on
        length_plot=cell.empty; 
        for i3=1:r3
            t=text_length(i3); 
            length_plot{i3}=num2str(t); 
        end 

        legend (length_plot, 'Location', 'EastOutside');  
        hold off 
        ylim([0 100]); 
        ylabel(['% of ' NAME ' with chain length']); 
        if T==6
            set(gca, 'XTickLabel', {'8', '16', '24', '32', '40', '48'}); 
        elseif T==3
            set (gca, 'XTickLabel', {'Control', 'Early', 'Late'}); 
        end 
        print (gcf, '-dpng', [label '_' NAME '_length.png']); 
        print (gcf, '-depsc2', [label '_' NAME '_length.eps']);
        close 

        [r3,c3]=size(plot_sat); 
        plot_sat=rot90(plot_sat); 
        plot_sat=flipud(plot_sat); 
        bar (1:c3, plot_sat, 'stack'); 
        hold on
        length_sat=cell.empty; 
        for i3=1:r3
            t=num2str(text_sat(i3)); 
            length_sat{i3}=t; 
        end 
        legend (length_sat, 'Location', 'EastOutside'); 
        ylim([0 100]); 
        ylabel(['% of ' NAME ' with saturation level']); 
        
        if T==6
            set(gca, 'XTickLabel', {'8', '16', '24', '32', '40', '48'}); 
        elseif T==3
            set (gca, 'XTickLabel', {'Control', 'Early', 'Late'}); 
        end 
        print (gcf, '-dpng', [label '_' NAME '_saturation.png']); 
        print (gcf, '-depsc2', [label '_' NAME '_saturation.eps']);
        close 
    else %% if its the control---------------------------------------
        for k2=1:r1
            for j=1:6
                time_mat(k2,j)=nanmean(levels(k2,1:3));
            end 
        end 
        sum_=nansum(time_mat);
        unique_sat=unique(saturation); 
        plot_sat=zeros(length(unique_sat),T); 
        text_sat=double.empty; 
        text_length=double.empty; 
        for j=1:length(unique_sat)
            u=unique_sat(j); 
            indx=find(saturation==u); %find the same saturation level
            lipid_new=time_mat(indx,:);
            [rl,cl]=size(lipid_new); 
            if rl>1
                lipid_sum=nansum(lipid_new); 
            else 
                lipid_sum=lipid_new; 
            end 
            for j2=1:6
                plot_sat(j,j2)=(lipid_sum(j2)/sum_(j2))*100;  
            end 
            text_sat(j)=u;
        end 
        unique_length=unique(chain_length); 
        plot_length=zeros(length(unique_length),T); 
        for j=1:length(unique_length) 
            u=unique_length(j); 
            indx=find(chain_length==u); %find the same saturation level
            lipid_new=time_mat(indx,:);
            [rl,c1]=size(lipid_new); 
            if rl>1
                lipid_sum=nansum(lipid_new); 
            else 
                lipid_sum=lipid_new; 
            end 
            for j2=1:6
                plot_length(j,j2)=(lipid_sum(j2)/sum_(j2))*100;  
            end 
            text_length(j)=u; 
        end 
        plot_combo=double.empty;
        count=1; 
        for j=1:length(unique_sat)
            us=unique_sat(j); 
            indxs=find(saturation==us); 
            for j2=1:length(unique_length)
                error=0; 
                ul=unique_length(j2); 
                indxl=find(chain_length==ul); 
                indx=intersect(indxl,indxs); 
                if numel(indx)==0 
                    error=1; 
                elseif numel(indx)>1
                    lipid_new=nansum(time_mat(indx,:));
                else 
                    lipid_new=time_mat(indx,:); 
                end 
                if error==0 
                    [rl,c1]=size(lipid_new); 
                    if rl>1
                        lipid_sum=nansum(lipid_new); 
                    else 
                        lipid_sum=lipid_new; 
                    end 
                    for j3=1:6
                        plot_combo(count,j3)=(lipid_sum(j3)/sum_(j3))*100;  
                    end 
                    text_combo(count,1)=us;
                    text_combo(count,2)=ul; 
                    count=count+1; 
                end
            end
        end 
        subplot(1,3,1);
       [r3,c3]=size(plot_combo);
        plot_combo=rot90(plot_combo); 
        plot_combo=flipud(plot_combo); 
        NAME=group{i}; 
        bar (1:2:c3*2, plot_combo, 'stack');
        hold on 
        combo_text=cell.empty; 
        for i3=1:r3
            t=num2str(text_combo(i3,1)); 
            v=num2str(text_combo(i3,2)); 
            x=[v ':' t]; 
            combo_text{i3}=x;  
        end 
        legend(combo_text, 'Location', 'EastOutside');    
        ylim([0 100]); 
        xlim([0 2]); 
        ylabel(['% of ' NAME ' with length:saturation']);    
        if T==6
            set(gca, 'XTickLabel', {'Control', '16', '24', '32', '40', '48'}); 
        elseif T==3
            set (gca, 'XTickLabel', {'Control', 'Early', 'Late'});  
        end 
        print (gcf, '-dpng', [label '_' NAME '_combo.png']); 
        print (gcf, '-depsc2', [label '_' NAME '_combo.eps']);
        close 

        subplot(1,3,1); 
        [r3,c3]=size(plot_length);
        plot_length=rot90(plot_length); 
        plot_length=flipud(plot_length);  
        bar (1:2:c3*2, plot_length, 'stack'); 
        hold on
        length_plot=cell.empty; 
        for i3=1:r3
            t=num2str(text_length(i3)); 
            length_plot{i3}=t; 
        end 

        legend (length_plot, 'Location', 'EastOutside');  
        hold off 
        ylim([0 100]); 
        ylabel(['% of ' NAME ' with chain length']); 
        if T==6
            set(gca, 'XTickLabel', {'Control', '16', '24', '32', '40', '48'}); 
        elseif T==3
            set (gca, 'XTickLabel', {'Control', 'Early', 'Late'});  
        end 
        xlim([0 2]);  
        print (gcf, '-dpng', [label '_' NAME '_length.png']); 
        print (gcf, '-depsc2', [label '_' NAME '_length.eps']);
        close 

        subplot(1,3,1); 
        [r3,c3]=size(plot_sat); 
        plot_sat=rot90(plot_sat); 
        plot_sat=flipud(plot_sat); 
        bar (1:2:c3*2, plot_sat, 'stack'); 
        length_sat=cell.empty; 
        for i3=1:r3
            t=num2str(text_sat(i3)); 
            length_sat{i3}=t; 
        end 
        legend (length_sat, 'Location', 'EastOutside'); 
        ylim([0 100]); 
        xlim([0 2]); 
        ylabel(['% of ' NAME ' with saturation level']); 
        if T==6
            set(gca, 'XTickLabel', {'Control', '16', '24', '32', '40', '48'}); 
        elseif T==3
            set (gca, 'XTickLabel', {'Control', 'Early', 'Late'});  
        end 
        print (gcf, '-dpng', [label '_' NAME '_saturation.png']); 
        print (gcf, '-depsc2', [label '_' NAME '_saturation.eps']);
        close 
    end
end 
end 


