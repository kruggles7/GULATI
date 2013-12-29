function [ anova_mat, ttest_mat, statdata1] = anova_ttest_final( mat, text_mat, time_points, samples )
%Uses normalized or nonnormalized lipidomic data and determines the
%significance across time points.  if the anova is significant it does a
%ttest to determine which time points are significant and uses the
%follwoing variables 
%   MAT is the data matrix with non_normalized values 
%THIS IS FOR 9 SAMPLES 6 TIME POINTS--CHANGE IF DIFFERENT

%OUTPUTS 
%   anova_mat : has the ANOVA p-values 
%   ttest_mat: has the significant pairs  

[row,col]= size(mat); 
test_num=row; 
anova_mat = zeros (test_num, 1); 
ttest_mat = cell (test_num, 3); 
statdata1=zeros(samples,time_points,row);  

for j=1:row 
    data=mat(j,:); 
    statdata=zeros(samples,6); 
    statdata(:,1)=data(1:samples); 
    statdata(:,2)=data(samples+1:2*samples); 
    statdata(:,3)=data((2*samples+1):3*samples); 
    statdata(:,4)=data((3*samples+1):4*samples); 
    statdata(:,5)=data((4*samples+1):5*samples); 
    statdata(:,6)=data((5*samples+1):6*samples);   
    if max(statdata)>0
        [p,a,s]=anova1(statdata);  
        %p=1; 
        anova_mat(j,1)=p; 
        ttest_mat {j,1}=text_mat(j,1); 
        ttest_mat{j,2}=p; 
        temp_ttest=[]; 
        %if p<0.05
            for n=1:time_points 
                group1=statdata(:,n); 
                for m=(n+1):time_points 
                    group2=statdata(:,m); 
                    [hval,pval]=ttest2(group1,group2); 
                    if pval<0.05
                        temp_ttest=[temp_ttest n m]; 
                    end 
                end 
            end 
            L=length(temp_ttest); 
            temp3=''; 
            for J=1:2:L
                temp1=num2str(temp_ttest(J)); 
                temp2=num2str(temp_ttest(J+1)); 
                temp3=[temp3 temp1 '/' temp2 ',']; 
            end 
            ttest_mat{j,3}=temp3; 
            if isempty(temp3)==1
                ttest_mat{j,3}='NA'; 
            end 
        %else 
         %  ttest_mat{j,3}='NA'; 
        %end 

        for m=1:time_points 
            A=[]; 
            for s=1:samples
                if isnan(statdata(s,m))==0
                    A=[A statdata(s,m)]; 
                else 
                    A=[A 0]; 
                end 
            end 
        end
        close all
    end 
    statdata1(:,:,j)=statdata; 
end 

end

