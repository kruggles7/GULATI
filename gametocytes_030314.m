
%used concentration values-----------------
cd ..
cd matrices 
load gametocytes_update
load time4
load time4_2
load species
cd ..
cd programs
g_mat(:, 1:45)=gametocyte(:, 1:45); 
g_mat(:, 46:54)=gametocyte(:, 49:57); 
[G_mat, mat_G_NaN, N]=filter_std(g_mat,9,6); 
l_malaria=G_mat(:,1:9); 
l_infected=G_mat(:,10:18); 
l_uninf=G_mat(:,19:27); 
e_malaria=G_mat(:,28:36); 
e_infected=G_mat(:,37:45); 
e_uninf=G_mat(:,46:54); 
late=G_mat(:,1:27); 
early=G_mat(:,28:54); 

late_NaN=mat_G_NaN(:,1:27); 
early_NaN=mat_G_NaN(:,28:54); 
[T4_mat, mat_T_NaN, n1]=filter_std(time4,3,1); 
[T4_mat2,mat_T_NaN_2, n2]=filter_std(time4_2,3,1); 

%[ttest, plot_mat]=Figure1_gametocytes( early, late, gametocyte_text,  species );
%[ttest_mat, plot_mat] = Figure1_gametocytes( early, late, g_names, species )

% [malaria_ttest, n3]= ttest_gametocyte( e_malaria, l_malaria, gametocyte_text, time4, time4_cell, time4_2, time4_2_cell); 
% [infected_ttest, n4]=ttest_gametocyte( e_infected, l_infected, gametocyte_text, T4_mat, time4_cell, time4_2, time4_2_cell); 
% [uninf_ttest, n5]=ttest_gametocyte(e_uninf, l_uninf, gametocyte_text, time4, time4_cell, time4_2, time4_2_cell); 
% [early_ttest]=ttest_gametocyte2(e_malaria, e_infected, e_uninf, gametocyte_text); 
% [late_ttest]=ttest_gametocyte2(l_malaria, l_infected, l_uninf, gametocyte_text); 
% 
% [malaria_ttest2]= ttest_gametocyte_trop1(e_malaria, l_malaria, gametocyte_text, T4_mat, time4_cell); 
% [infected_ttest2]=ttest_gametocyte_trop1(e_infected, l_infected, gametocyte_text, T4_mat, time4_cell); 
% [uninf_ttest2]=ttest_gametocyte_trop1(e_uninf, l_uninf, gametocyte_text, T4_mat, time4_cell); 

% [r,c]=size(early); 
% %heatmap 
% [ plot_mat, stat_mat, stat_mat2] = gametocyte_heatmap( e_malaria(14:r,:), l_malaria(14:r,:), gametocyte_text(14:r,:), time4, time4_2, time4_cell, time4_2_cell, species); 
% %[stat_mat, plot_mat]=Figure1_gametocytes( early, late, gametocyte_text, species ); 
% 
% H=HeatMap(HM_names2, 'Colormap', cmap, 'Symmetric','false'); 

% [plot_mat, HM_names2]= gametocyte_heatmap_step2 ( e_malaria(14:r,:), l_malaria(14:r,:), gametocyte_text(14:r,:),time4, time4_cell,time4_2, time4_2_cell ); 
%  H=HeatMap(HM_names2, 'Colormap', cmap, 'Symmetric','false'); 
% plot(H);
%  cd ..
%  cd results
%  cd heatmaps
%  print (gcf, '-dpng', 'gametocyte_names.png')
%  print (gcf, '-depsc2', ['gametocyte_names.eps']); 
%  %saveas (gcf, 'parasite_names.fig'); 
%  close all
%  cd ..
%  cd .. 
%  cd programs

% both=cat(2,early,late); 
% [r,c]=size(both); 
% both=both(14:r, :); 
% gam_text=gametocyte_text(14:r,:); 
% [anova_sat_G, anova_chain_G, anova_combo_G, plot_combo_S, combo_text_S]= sat_length_gams(both, gam_text, species, 'gams', 2); 
