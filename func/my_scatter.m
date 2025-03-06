function my_scatter(x, y, title_label, x_label, y_label)

    scatter(x,y, 'o', 'filled');
    [rho_s, p_s] = corr(x,y, 'Type', 'S', 'rows', 'complete');
    title(sprintf('%s: r = %.2f p = %f', title_label, rho_s, p_s/2), 'FontWeight', 'normal'); % one-tailed test
    xlabel(x_label);
    ylabel(y_label);
    set(gca, 'FontSize', 14);
    
end