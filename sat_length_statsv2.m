function [anova_sat_RC,anova_sat_PR,anova_sat_PC, anova_length_RC,anova_length_PR, anova_length_PC, anova_combo_RC,anova_combo_PR, anova_combo_PC] = sat_length_statsv2( mat_control, matp, matR, lipids, species )
%UNTITLED Summary of this function goes here
% %   Detailed explanation goes here
% % 
%  mat_control=uninf;  
%  matp=DD2f_p_; 
%  matR=DD2f_R_; 
%  lipids=DD2f_names; 
%  species=species; 
 
[r,c]=size(lipids); 
group={'PL', 'SL', 'NL'}; 
match=cell2mat(species(:,3)); 
anova_sat=cell.empty;
anova_chain=cell.empty; 
anova_combo=cell.empty; 
a_count=1; 
l_count=1; 
co_count=1; 

match_=match; 
match_(match==4)=3; 
for i=1:3
    indx=find(match_==i); 
    L=species(indx,2); 
    saturation=double.empty; 
    chain_length=double.empty;
    sum_=double.empty; 
    time_mat=double.empty; 
    col_label=double.empty; 
    levels_control=double.empty; 
    levels_p=double.empty; 
    levels_R=double.empty; 
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
                    levels_control(n_,:)=mat_control(j,:); 
                    levels_p(n_,:)=matp(j,:); 
                    levels_R(n_,:)=matR(j,:); 
                    col_label{n_,1}=text1;  
                    n_=n_+1;
                end  
            end 
        end
    end  
    %make an average with time for each member of the group
    %control--------------------------------------------------------------
    [r1,c1]=size(levels_control); 
    levels_sum=sum(levels_control); 
    levels_percent_c=zeros(r1,c1); 
    for k2=1:r1
        for j2=1:c1
            levels_percent_c(k2,j2)=(levels_control(k2,j2)/levels_sum(j2))*100; 
        end 
    end 
    %parasite-------------------------------------------------------------
    [r1,c1]=size(levels_p); 
    levels_sum=sum(levels_p); 
    levels_percent_p=zeros(r1,c1); 
    for k2=1:r1
        for j2=1:c1
            levels_percent_p(k2,j2)=(levels_p(k2,j2)/levels_sum(j2))*100; 
        end 
    end 
    %RBC-------------------------------------------------------------------
    [r1,c1]=size(levels_R); 
    levels_sum=sum(levels_R); 
    levels_percent_R=zeros(r1,c1); 
    for k2=1:r1
        for j2=1:c1
            levels_percent_R(k2,j2)=(levels_R(k2,j2)/levels_sum(j2))*100; 
        end 
    end 
    unique_sat=unique(saturation); 
    stat_mat_RC=double.empty; 
    stat_mat_PC=double.empty; 
    stat_mat_PR=double.empty; 
    %for each saturation---------------------------------------------------
    for j=1:length(unique_sat) 
        u=unique_sat(j); 
        anova_sat_RC{a_count,1}=group{i};
        anova_sat_RC{a_count,2}=u; 
        anova_sat_PR{a_count,1}=group{i};
        anova_sat_PR{a_count,2}=u; 
        anova_sat_PC{a_count,1}=group{i};
        anova_sat_PC{a_count,2}=u; 
        indx=find(saturation==u);
        stat_sat_c=levels_percent_c(indx,:); 
        stat_sat_R=levels_percent_R(indx,:); 
        stat_sat_p=levels_percent_p(indx,:); 
        [r3,c3]=size(stat_sat_c); 
        if r3>1
            stat_sat_c=sum(stat_sat_c); 
        end  
        [r3,c3]=size(stat_sat_R); 
        if r3>1
            stat_sat_R=sum(stat_sat_R); 
        end
        [r3,c3]=size(stat_sat_p); 
        if r3>1
            stat_sat_p=sum(stat_sat_p); 
        end 
        d=1; 
        %for each timepoint-----------------------------------------------
        for k=1:6 %for each timepoint 
            temp=stat_sat_p(d:d+8); 
            group_p=temp;
            temp=stat_sat_R(d:d+8); 
            group_R=temp ; 
            d=d+9; 
            temp=stat_sat_c(:); 
            group_c=temp;  
            [hval, pval]=ttest2(group_p,group_R, 0.05,'both', 'unequal'); 
            anova_sat_PR{a_count,k+2}=pval; 
            [hval, pval]=ttest2(group_p,group_c, 0.05,'both', 'unequal');
            anova_sat_PC{a_count,k+2}=pval; 
            [hval, pval]=ttest2(group_c,group_R, 0.05,'both','unequal');
            anova_sat_RC{a_count,k+2}=pval; 
        end 
        a_count=a_count+1; 
    end 
    %length---------------------------------------------------------------
     unique_length=unique(chain_length);
    for j=1:length(unique_length) 
        u=unique_length(j); 
        anova_length_RC{l_count,1}=group{i};
        anova_length_RC{l_count,2}=u; 
        anova_length_PR{l_count,1}=group{i};
        anova_length_PR{l_count,2}=u; 
        anova_length_PC{l_count,1}=group{i};
        anova_length_PC{l_count,2}=u; 
        indx=find(chain_length==u);
        stat_length_c=levels_percent_c(indx,:); 
        stat_length_R=levels_percent_R(indx,:); 
        stat_length_p=levels_percent_p(indx,:); 
        [r3,c3]=size(stat_length_c); 
        if r3>1
            stat_length_c=sum(stat_length_c); 
        end  
        [r3,c3]=size(stat_length_R); 
        if r3>1
            stat_length_R=sum(stat_length_R); 
        end 
        [r3,c3]=size(stat_length_p); 
        if r3>1
            stat_length_p=sum(stat_length_p); 
        end 
        d=1; 
        %for each timepoint-----------------------------------------------
        for k=1:6 %for each timepoint 
            temp=stat_length_p(d:d+8); 
            group_p=temp;
            temp=stat_length_R(d:d+8); 
            group_R=temp ; 
            d=d+9; 
            temp=stat_length_c(:); 
            group_c=temp;  
            [hval, pval]=ttest2(group_p,group_R, 0.05,'both', 'unequal');
            anova_length_PR{l_count,k+2}=pval; 
            [hval, pval]=ttest2(group_p,group_c, 0.05,'both', 'unequal'); 
            anova_length_PC{l_count,k+2}=pval; 
            [hval, pval]=ttest2(group_c,group_R, 0.05,'both', 'unequal'); 
            anova_length_RC{l_count,k+2}=pval; 
        end 
        l_count=l_count+1; 
    end  
    %combo---------------------------------------------------------------
    for j=1:length(unique_sat)
        us=unique_sat(j); 
        indxs=find(saturation==us); 
        for j2=1:length(unique_length)
            error=0; 
            ul=unique_length(j2); 
            indxl=find(chain_length==ul); 
            indx=intersect(indxl,indxs); 
            if numel(indx)>0
                anova_combo_RC{co_count,1}=group{i};
                anova_combo_RC{co_count,2}=us; 
                anova_combo_RC{co_count,3}=ul; 
                anova_combo_PR{co_count,1}=group{i};
                anova_combo_PR{co_count,2}=us; 
                anova_combo_PR{co_count,3}=ul; 
                anova_combo_PC{co_count,1}=group{i};
                anova_combo_PC{co_count,2}=us; 
                anova_combo_PC{co_count,3}=ul;
                stat_length_c=levels_percent_c(indx,:); 
                stat_length_R=levels_percent_R(indx,:); 
                stat_length_p=levels_percent_p(indx,:); 
                [r3,c3]=size(stat_length_c); 
                if r3>1
                    stat_length_c=sum(stat_length_c); 
                end  
                [r3,c3]=size(stat_length_R); 
                if r3>1
                    stat_length_R=sum(stat_length_R); 
                end 
                [r3,c3]=size(stat_length_p); 
                if r3>1
                    stat_length_p=sum(stat_length_p); 
                end
                d=1; 
                %for each timepoint-----------------------------------------------
                for k=1:6 %for each timepoint 
                    temp=stat_length_p(d:d+8); 
                    group_p=temp;
                    temp=stat_length_R(d:d+8); 
                    group_R=temp ; 
                    d=d+9; 
                    temp=stat_length_c(:); 
                    group_c=temp;  
                    [hval, pval]=ttest2(group_p,group_R, 0.05,'both', 'unequal'); 
                    anova_combo_PR{co_count,k+3}=pval; 
                    [hval, pval]=ttest2(group_p,group_c, 0.05,'both', 'unequal'); 
                    anova_combo_PC{co_count,k+3}=pval; 
                    [hval, pval]=ttest2(group_c,group_R, 0.05,'both', 'unequal');
                    anova_combo_RC{co_count,k+3}=pval; 
                end 
                co_count=co_count+1; 
            end 
        end 
    end 
end 
end 




