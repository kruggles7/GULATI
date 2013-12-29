function [ ] = plot_combo( combo_p, combo_R,combo_u, combo_text )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[r,c]=size(combo_p); 
for i=1%:c
    Tp=combo_p(:,i); 
    TR=combo_R(:,i); 
    Tu=combo_u(:,i); 
    name=combo_text{i}; 
    x=1:6; 
    plot(x, Tp,'-sb', 'LineWidth', 1.5, 'MarkerFaceColor', 'b'); 
    hold on
    plot(x, TR, '-ok', 'LineWidth', 1.5, 'MarkerFaceColor', 'k'); 
    plot(x, Tu, '--r', 'LineWidth',2); 
    hold off
    title(name); 
    set(gca, 'XTick', [0 1 2 3 4 5 6]); 
    set (gca, 'XTickLabel', [8 16 24 32 40 48]); 
    xlabel('Time (hours)'); 
    ylabel('% Total Species'); 
end 

end

