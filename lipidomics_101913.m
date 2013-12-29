% 
%create matrices from cell arrays DD2-------------
load DD2_names.mat
load PFS_DD2.mat 
load species.mat

DD2f=DD2fraction(4:307,2:109); 
DD2f_parasite=cell2mat(DD2fraction(4:307,2:55)); 
DD2f_RBC=cell2mat(DD2fraction(4:307,56:109)); 
DD2f_p=reorder(DD2f_parasite);
DD2f_R=reorder(DD2f_RBC); 
[DD2f_p_, mat_p_NaN]=filter_std(DD2f_p,9,6); 
[DD2f_R_, mat_R_NaN]=filter_std(DD2f_R,9,6); 

%uninfected RBC
uninfected=cell2mat(DD2fraction(4:307, 113:115));
[uninf, uinf_NaN]=filter_std(uninfected,3,1); 

%FIGURE 1- lipid species bar charts
%[mol_parasite, mol_rbc, mol_uninfected, plot_mat]=Figure1( DD2f_p_, DD2f_R_, uninf, DD2f_names, species );
%[ anova_parasite, ttest_parasite, statdatap] = anova_ttest_final( DD2f_p_, DD2f_names, 6, 9 ); 
%[ anova_RBC, ttest_RBC, statdataR] = anova_ttest_final( DD2f_R_, DD2f_names, 6, 9 ); 
%[ pvalue_p_R, pvalue_p_u, pvalue_R_u ] = ttest_within_species( DD2f_p_, DD2f_R_, uninf, DD2f_names, 6); 

%[anova_parasite_mol, ttest_parasite_mol, statdatatp_mol]=anova_ttest_final(mol_parasite, species(:,2), 6,9);
%[anova_RBC_mol, ttest_RBC_mol, statdataR_mol]=anova_ttest_final(mol_rbc, species(:,2), 6,9); 
%[pvalue_p_R_mol, pvalue_p_u_mol, pvalue_R_u_mol]= ttest_within_species(mol_parasite, mol_rbc, mol_uninfected, species(:,2), 6); 
 
% %FIGURE 2 - heatmap
% [ plot_mat_R, bar_mat_R, percent_changed_R] = lipid_heatmap ( DD2f_R_, DD2f_names, 'RBC', species ); 
% [ plot_mat_p, bar_mat_p, percent_changed_p ] = lipid_heatmap ( DD2f_p_, DD2f_names, 'parasite',species ); 
% 
% 
%length and saturation bars 
 %[anova_sat_p, anova_chain_p, anova_combo_p, plot_combo_p, combo_text_p]= sat_length(DD2f_p_, DD2f_names, species, 'parasite', 6); 
 %[anova_sat_R, anova_chain_R, anova_combo_R, plot_combo_R, combo_text_R]= sat_length(DD2f_R_, DD2f_names, species, 'RBC', 6);
 %[anova_sat_u, anova_chain_u, anova_combo_u, plot_combo_u, combo_text_u]=sat_length(uninf, DD2f_names, species, 'uninfected',1); 

%plot combo info
%plot_combo(plot_combo_p, plot_combo_R, plot_combo_u, combo_text_p); 

%supernatant--------------------------------------------------
 load super.mat 
 
 super_mat=cell2mat(Supp_species(4:292,2:55)); 
 super_names=Supp_species(4:292,1); 
 [super_, super_NaN]=filter_std(super_mat,9, 6); 
 [anova_super, ttest_super, statdata_super]=anova_ttest_final(super_, super_names, 6, 9); 
 %[plot_super, bar_mat_super, percent_super]=lipid_heatmap(super_,super_names,'super', species); 
% 
