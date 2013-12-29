% 
%create matrices from cell arrays DD2-------------
load DD2_names.mat
load PFS_DD2.mat 
DD2f=DD2fraction(4:307,2:109); 
%DD2f_names=DD2fraction(4:307,1); 
DD2f_parasite=cell2mat(DD2fraction(4:307,2:55)); 
DD2f_RBC=cell2mat(DD2fraction(4:307,56:109)); 
DD2f_p=reorder(DD2f_parasite);
DD2f_R=reorder(DD2f_RBC); 
[DD2f_p_, mat_p_NaN]=filter_std(DD2f_p,9); 
[DD2f_R_, mat_R_NaN]=filter_std(DD2f_R,9); 

[ plot_mat_R ] = lipid_heatmap ( DD2f_R_, DD2f_names, 'RBC' ); 
[ plot_mat_p ] = lipid_heatmap ( DD2f_p_, DD2f_names, 'parasite' ); 


%supernatant--------------------------------------------------
load super.mat 

super_mat=cell2mat(Supp_species(4:292,2:55)); 
super_names=Supp_species(4:292,1); 
[super_, super_NaN]=filter_std(super_mat,9); 
[plot_super]=lipid_heatmap(super_,super_names,'super'); 

