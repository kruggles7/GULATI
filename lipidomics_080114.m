
%create matrices from cell arrays DD2-------------
cd ..
cd matrices
load DD2_names.mat
load PFS_DD2.mat 
load species.mat
cd ..
cd programs

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

%  %FIGURE 1- lipid species bar charts
[mol_parasite1,  mol_uninfecte1, plot_mat1, mol_names1]=Figure1_noinfected( DD2f_p_, uninf, DD2f_names, species );
[mol_parasite, mol_rbc, mol_uninfected, plot_mat, mol_names, mol_parasite_all, mol_rbc_all, mol_uninf_all ]=Figure1_v2 (DD2f_p_, DD2f_R_, uninf, DD2f_names, species ); 

%FIGURE 1- lipid species bar charts
[ anova_parasite, ttest_parasite, statdatap] = anova_ttest_final( DD2f_p_, DD2f_names, 6, 9 ); 
[ anova_RBC, ttest_RBC, statdataR] = anova_ttest_final( DD2f_R_, DD2f_names, 6, 9 ); 
[ pvalue_p_R, pvalue_p_u, pvalue_R_u ] = ttest_within_species( DD2f_p_, DD2f_R_, uninf, DD2f_names, 6); 

[anova_parasite_mol, ttest_parasite_mol, statdatatp_mol]=anova_ttest_final(mol_parasite, mol_names1, 6,9);
[anova_RBC_mol, ttest_RBC_mol, statdataR_mol]=anova_ttest_final(mol_rbc, mol_names1, 6,9); 
[pvalue_p_R_mol, pvalue_p_u_mol, pvalue_R_u_mol]= ttest_within_species(mol_parasite, mol_rbc, mol_uninfected, mol_names1, 6); 
 

%  cd ..
%  cd results
%  mkdir ('heatmaps')
%  cd ..
%  cd programs
% %FIGURE 2 - heatmap
%  [ plot_mat_R, bar_mat_R, percent_changed_R, stats_R, names_R] = lipid_heatmap ( DD2f_R_, DD2f_names, 'RBC', species); 
%  cmap='jet';
%  R=flipud(names_R); 
%  H=HeatMap(R, 'Colormap', cmap, 'Symmetric','false' );
%  plot(H); 
%  cd ..
%  cd results
%  cd heatmaps
%  print (gcf, '-dpng', 'RBC_names.png'); 
%  print (gcf, '-depsc2', ['RBC_names.eps']);
%  %saveas (gcf, 'RBC_names.fig'); 
%  cd ..
%  cd .. 
%  cd programs
%  close all
 % [ plot_mat_p, bar_mat_p, percent_changed_p, stats_p, names_p ] = lipid_heatmap ( DD2f_p_, DD2f_names, 'parasite', species); 
%  P=flipud(names_p); 
%   H=HeatMap(P, 'Colormap', cmap, 'Symmetric','false'); 
%  plot(H);
%  cd ..
%  cd results
%  cd heatmaps
%  print (gcf, '-dpng', 'parasite_names.png')
%  print (gcf, '-depsc2', [ 'parasite_names.eps']); 
%  %saveas (gcf, 'parasite_names.fig'); 
%  close all
%  cd ..
%  cd .. 
%  cd programs

 
%length and saturation bars 
%  [anova_sat_p, anova_chain_p, anova_combo_p]= sat_length_simple(DD2f_p_, DD2f_names, species, 'parasite', 6, 0); 
%  cd ..
%  cd ..
%  cd programs
%  [anova_sat_R, anova_chain_R, anova_combo_R]= sat_length_simple(DD2f_R_, DD2f_names, species, 'RBC', 6, 0);
%  cd ..
%  cd ..
%  cd programs
%  [anova_sat_u, anova_chain_u, anova_combo_u]=sat_length_simple(uninf, DD2f_names, species, 'uninfected',1, 0); 
%  cd ..
%  cd ..
%  cd programs
 
 %[anova_sat_RC,anova_sat_PR,anova_sat_PC, anova_length_RC,anova_length_PR, anova_length_PC, anova_combo_RC,anova_combo_PR, anova_combo_PC] = sat_length_stats( uninf, DD2f_p_, DD2f_R_, DD2f_names, species ); 


% %supernatant--------------------------------------------------
% cd ..
% cd matrices
% load super.mat 
% load Supp_control.mat
% load super_names; 
% cd ..
% cd programs
%  
%  super_mat=cell2mat(Supp_species(4:292,2:55));
%  [super_, super_NaN]=filter_std(super_mat,9, 6);
%  [super_cont, super_cont_NaN]=filter_std(control, 9, 2); 

%  [mol_super,mol_control, mol_uiRBC, plot_mat, mol_names]=Figure1_super( super_, control, super_names, species );
%  [anova_super_mol, ttest_super_mol, statdatatS_mol]=anova_ttest_final(mol_super, mol_names, 6,9);
%  [ pvalue_s_c, pvalue_s_u] = ttest_within_species_super(mol_super, mol_control, mol_uiRBC, mol_names, 6); 

 
% [anova_super, ttest_super, statdata_super]=anova_ttest_final(super_, super_names, 6, 9); 

% [plot_super, bar_mat_super, percent_super, stats_super, names_S]=lipid_heatmap(super_,super_names,'super', species); 
%  cmap='jet';
%  S=flipud(names_S); 
%  H=HeatMap(S,'Colormap', cmap, 'Symmetric','false'); 
% plot(H);
%  cd ..
%  cd results
%  cd heatmaps
%  print (gcf, '-dpng', 'super_names.png')
%  print (gcf, '-depsc2','super_names.eps'); 
%  %saveas (gcf, 'super_names.fig')
%  cd ..
%  cd .. 
%  cd programs
%  close all
%  
% [anova_sat_S, anova_chain_S, anova_combo_S]= sat_length_simple(super_, super_names, species, 'super', 6, 0); 

