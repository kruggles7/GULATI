function [anova_sat, anova_chain, anova_combo ] = sat_length_simple( mat, lipids, species, label, T, stat_yes )
%UNTITLED Summary of this function goes here
% %   Detailed explanation goes here
% 
mat=DD2f_p_; 
lipids=DD2f_names; 
label='parasite'; 
T=6; 
stat_yes=0

[r,c]=size(lipids); 
group={'PL', 'SL', 'NL', 'DAG', 'APG'}; 
match_=cell2mat(species(:,3)); 
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
for i=1%:5
    indx=find(match_==i); 
    L=species(indx,2); 
    saturation=double.empty; 
    chain_length=double.empty;
    sum_=double.empty; 
    time_mat=double.empty; 
    col_label=double.empty; 
    levels=double.empty; 
    n_=1; 
    NAME=group{i}; 
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
                        rL=P-1; 
                    else 
                        rL=length(text_);
                    end 
                    for p_=q1-1:q1
                        c_=text_(p_);  
                        chain=[chain c_] ;
                    end
                    chain = str2num(chain); 
                    %chain_length(n_,1)=chain; 
                    
                    for p_=q1+2:rL
                        c_=text_(p_) ;
                        sat_=[sat_ c_]; 
                    end 
                    sat_=str2num(sat_);
               %     saturation(n_,1)=sat_; 
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
    levels_sum=sum(levels); 
    levels_percent=zeros(r1,c1); 
    for k2=1:r1
        for j2=1:c1
            levels_percent(k2,j2)=(levels(k2,j2)/levels_sum(j2))*100; 
        end 
    end 
    if T>1 %% ALL EXCEPT FOR CONTROL 
        for k2=1:r1
            d=1; 
            for j=1:T 
                time_mat(k2,j)=mean(levels(k2,d:d+8)); 
                d=d+9; 
            end
        end
        sum_=sum(time_mat);
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
                lipid_sum=sum(lipid_new); 
            else 
                lipid_sum=lipid_new; 
            end 
            for j2=1:T
                plot_sat(j,j2)=(lipid_sum(j2)/sum_(j2))*100;  
            end 
            text_sat(j)=u; 
            stat_sat=levels_percent(indx,:); 
            stat_mat=double.empty; 
            [r3,c3]=size(stat_sat); 
            if r3>1
                stat_sat=sum(stat_sat); 
            end  
            d=1; 
            for j2=1:T 
                temp=stat_sat(d:d+8); 
                stat_mat(:,j2)=temp; 
                d=d+9; 
            end
            if stat_yes==1
                p_sat=anova1(stat_mat); 
                close all
            else 
                p_sat='NA'; 
            end 
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
                lipid_sum=sum(lipid_new); 
            else 
                lipid_sum=lipid_new; 
            end 
            for j2=1:T
                plot_length(j,j2)=(lipid_sum(j2)/sum_(j2))*100;  
            end 
            text_length(j)=u; 
            stat_chain=levels_percent(indx,:); 
            stat_mat=double.empty; 
            [r3,c3]=size(stat_chain); 
            if r3>1
                stat_chain=sum(stat_chain); 
            end  
            d=1; 
            for j2=1:T 
                temp=stat_chain(d:d+8); 
                stat_mat(:,j2)=temp; 
                d=d+9; 
            end
            if stat_yes==1
                p_chain=anova1(stat_mat); 
                close all
            else 
                p_chain='NA'; 
            end 
            
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
                    lipid_new=sum(time_mat(indx,:));
                else 
                    lipid_new=time_mat(indx,:); 
                end 
                if error==0 
                    [rl,c1]=size(lipid_new); 
                    if rl>1
                        lipid_sum=sum(lipid_new); 
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
                    stat_mat=double.empty; 
                    [r3,c3]=size(stat_combo); 
                    if r3>1
                        stat_combo=sum(stat_combo); 
                    end  
                    d=1; 
                    for j2=1:T 
                        temp=stat_combo(d:d+8); 
                        stat_mat(:,j2)=temp; 
                        d=d+9; 
                    end
                    if stat_yes==1
                        p_chain=anova1(stat_mat); 
                        close all
                    else 
                        p_chain='NA'; 
                    end 
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
   
        
        %plot length
        [r3,c3]=size(plot_length); 
        plot_length_=zeros(c3,4);  
        for S=1:r3
            if text_length(S)<=12
                plot_length_(:,1)=plot_length_(:,1)+ plot_length(S,:)'; 
            elseif text_length(S)>12 && text_length(S)<=24
                plot_length_(:,2)=plot_length_(:,2)+ plot_length(S,:)'; 
            elseif text_length(S)<42
                 plot_length_(:,3)=plot_length_(:,3)+ plot_length(S,:)'; 
            else 
                plot_length_(:,4)=plot_length_(:,4)+plot_length(S,:)'; 
            end 
        end 
        
        f1=subplot(6,1,1);
        bar(1:4, plot_length_(1,:)); 
        set(gca, 'XTickLabel', ''); 
        ylim([0 100]); 
        ylabel('8'); 
        title(['% of ' NAME  ' with chain length']); 
        
        f2=subplot(6,1,2);
        bar(1:4,plot_length_(2,:)); 
        set(gca, 'XTickLabel', ''); 
        ylim([0 100]); 
        ylabel('16'); 
        linkaxes([f1 f2],'x'); %make y axis the same
        pos1=get(f1,'Position'); %find the current position [x,y,width,height]
        pos2=get(f2,'Position'); %find the current position [x,y,width,height]
        pos2(2)=pos1(2) - pos2(4); %move the second so it touches the first 
        set (f2,'Position',pos2); 
        
        f3=subplot(6,1,3); 
        bar(1:4,plot_length_(3,:)); 
        ylim([0 100]); 
        set(gca, 'XTickLabel', ''); 
        ylabel('24'); 
        linkaxes([f2 f3],'x'); %make y axis the same
        pos1=get(f2,'Position'); %find the current position [x,y,width,height]
        pos2=get(f3,'Position'); %find the current position [x,y,width,height]
        pos2(2)=pos1(2) - pos2(4); %move the second so it touches the first 
        set (f3,'Position',pos2); 
        
        f4=subplot(6,1,4); 
        bar(1:4,plot_length_(4,:)); 
        ylim([0 100]); 
        set(gca, 'XTickLabel', ''); 
        ylabel('32'); 
        linkaxes([f3 f4],'x'); %make y axis the same
        pos1=get(f3,'Position'); %find the current position [x,y,width,height]
        pos2=get(f4,'Position'); %find the current position [x,y,width,height]
        pos2(2)=pos1(2) - pos2(4); %move the second so it touches the first 
        set (f4,'Position',pos2); 
        
        f5=subplot(6,1,5); 
        bar(1:4,plot_length_(5,:)); 
        ylim([0 100]); 
        set(gca, 'XTickLabel', '');  
        ylabel('40'); 
        linkaxes([f4 f5],'x'); %make y axis the same
        pos1=get(f4,'Position'); %find the current position [x,y,width,height]
        pos2=get(f5,'Position'); %find the current position [x,y,width,height]
        pos2(2)=pos1(2) - pos2(4); %move the second so it touches the first 
        set (f5,'Position',pos2); 
        
        f6=subplot(6,1,6); 
        bar(1:4,plot_length_(6,:)); 
        ylim([0 100]); 
        set(gca, 'XTickLabel', {'Short', 'Medium', 'Long', 'Very Long'}); 
        ylabel('48'); 
        linkaxes([f5 f6],'x'); %make y axis the same
        pos1=get(f5,'Position'); %find the current position [x,y,width,height]
        pos2=get(f6,'Position'); %find the current position [x,y,width,height]
        pos2(2)=pos1(2) - pos2(4); %move the second so it touches the first 
        set (f6,'Position',pos2);
        
        
        print (gcf, '-dpng', [label '_' NAME '_length_simple.png']); 
        print (gcf, '-depsc2', [label '_' NAME '_length_simple.eps']);
        close 

        [r3,c3]=size(plot_sat); 
        plot_sat_=zeros(c3,3);  
        for S=1:r3
            if text_sat(S)==0
                plot_sat_(:,1)=plot_sat_(:,1)+ plot_sat(S,:)'; 
            elseif text_sat(S)==1
                plot_sat_(:,2)=plot_sat_(:,2)+ plot_sat(S,:)'; 
            elseif text_sat(S)>1
                 plot_sat_(:,3)=plot_sat_(:,3)+ plot_sat(S,:)'; 
            else 
                text_sat(S)
            end 
        end 
        
        f1=subplot(6,1,1);
        bar(1:3, plot_sat_(1,:)); 
        set(gca, 'XTickLabel', ''); 
        ylim([0 100]); 
        ylabel('8'); 
        title(['% of ' NAME ' with saturation level']); 
        
        f2=subplot(6,1,2);
        bar(1:3,plot_sat_(2,:)); 
        set(gca, 'XTickLabel', ''); 
        ylim([0 100]); 
        ylabel('16'); 
        linkaxes([f1 f2],'x'); %make y axis the same
        pos1=get(f1,'Position'); %find the current position [x,y,width,height]
        pos2=get(f2,'Position'); %find the current position [x,y,width,height]
        pos2(2)=pos1(2) - pos2(4); %move the second so it touches the first 
        set (f2,'Position',pos2); 
        
        f3=subplot(6,1,3); 
        bar(1:3,plot_sat_(3,:)); 
        ylim([0 100]); 
        set(gca, 'XTickLabel', ''); 
        ylabel('24'); 
        linkaxes([f2 f3],'x'); %make y axis the same
        pos1=get(f2,'Position'); %find the current position [x,y,width,height]
        pos2=get(f3,'Position'); %find the current position [x,y,width,height]
        pos2(2)=pos1(2) - pos2(4); %move the second so it touches the first 
        set (f3,'Position',pos2); 
        
        f4=subplot(6,1,4); 
        bar(1:3,plot_sat_(4,:)); 
        ylim([0 100]); 
        set(gca, 'XTickLabel', ''); 
        ylabel('32'); 
        linkaxes([f3 f4],'x'); %make y axis the same
        pos1=get(f3,'Position'); %find the current position [x,y,width,height]
        pos2=get(f4,'Position'); %find the current position [x,y,width,height]
        pos2(2)=pos1(2) - pos2(4); %move the second so it touches the first 
        set (f4,'Position',pos2); 
        
        f5=subplot(6,1,5); 
        bar(1:3,plot_sat_(5,:)); 
        ylim([0 100]); 
        set(gca, 'XTickLabel', '');  
        ylabel('40'); 
        linkaxes([f4 f5],'x'); %make y axis the same
        pos1=get(f4,'Position'); %find the current position [x,y,width,height]
        pos2=get(f5,'Position'); %find the current position [x,y,width,height]
        pos2(2)=pos1(2) - pos2(4); %move the second so it touches the first 
        set (f5,'Position',pos2); 
        
        f6=subplot(6,1,6); 
        bar(1:3,plot_sat_(6,:)); 
        ylim([0 100]); 
        set(gca, 'XTickLabel', {'SAT', 'MUFA', 'PUFA'}); 
        ylabel('48'); 
        linkaxes([f5 f6],'x'); %make y axis the same
        pos1=get(f5,'Position'); %find the current position [x,y,width,height]
        pos2=get(f6,'Position'); %find the current position [x,y,width,height]
        pos2(2)=pos1(2) - pos2(4); %move the second so it touches the first 
        set (f6,'Position',pos2); 
        print (gcf, '-dpng', [label '_' NAME '_saturation_simple.png']); 
        print (gcf, '-depsc2', [label '_' NAME '_saturation_simple.eps']);
        close 
    else %% if its the control---------------------------------------
        for k2=1:r1
            for j=1:6
                time_mat(k2,j)=mean(levels(k2,1:3));
            end 
        end 
        sum_=sum(time_mat);
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
                lipid_sum=sum(lipid_new); 
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
                lipid_sum=sum(lipid_new); 
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
                    lipid_new=sum(time_mat(indx,:));
                else 
                    lipid_new=time_mat(indx,:); 
                end 
                if error==0 
                    [rl,c1]=size(lipid_new); 
                    if rl>1
                        lipid_sum=sum(lipid_new); 
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
        %plot length
        [r3,c3]=size(plot_length); 
        plot_length_=zeros(c3,4);  
        for S=1:r3
            if text_length(S)<=12
                plot_length_(:,1)=plot_length_(:,1)+ plot_length(S,:)'; 
            elseif text_length(S)>12 && text_length(S)<=24
                plot_length_(:,2)=plot_length_(:,2)+ plot_length(S,:)'; 
            elseif text_length(S)<42
                 plot_length_(:,3)=plot_length_(:,3)+ plot_length(S,:)'; 
            else 
                plot_length_(:,4)=plot_length_(:,4)+plot_length(S,:)'; 
            end 
        end 
        
        f1=subplot(6,1,1);
        bar(1:4, plot_length_(1,:)); 
        set(gca, 'XTickLabel', {'Short', 'Medium', 'Long', 'Very Long'}); 
        ylim([0 100]); 
        ylabel('Uninfected'); 
        title(['% of ' NAME  ' with chain length']); 
        print (gcf, '-dpng', [label '_' NAME '_length.png']); 
        print (gcf, '-depsc2', [label '_' NAME '_length.eps']);
        close 

        [r3,c3]=size(plot_sat); 
        plot_sat_=zeros(c3,3);  
        for S=1:r3
            if text_sat(S)==0
                plot_sat_(:,1)=plot_sat_(:,1)+ plot_sat(S,:)'; 
            elseif text_sat(S)==1
                plot_sat_(:,2)=plot_sat_(:,2)+ plot_sat(S,:)'; 
            elseif text_sat(S)>1
                 plot_sat_(:,3)=plot_sat_(:,3)+ plot_sat(S,:)'; 
            else 
                text_sat(S)
            end 
        end 
        
        f1=subplot(6,1,1);
        bar(1:3, plot_sat_(1,:)); 
        set(gca, 'XTickLabel', {'SAT', 'MUFA', 'PUFA'}); 
        ylim([0 100]); 
        ylabel('Uninfected'); 
        title(['% of ' NAME ' with saturation level']); 
         
        print (gcf, '-dpng', [label '_' NAME '_saturation.png']); 
        print (gcf, '-depsc2', [label '_' NAME '_saturation.eps']);
        close 
    end
end 
end 


