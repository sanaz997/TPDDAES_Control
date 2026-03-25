function eig_plot_map(Map)
plot(real(Map(1,:)),imag(Map(1,:)), '-.','Color',[.3 .3 .3], 'MarkerSize', 20, 'LineWidth', 3)
hold on;
plot(real(Map(2,:)),imag(Map(2,:)), '-.','Color',[.3 .3 .3], 'MarkerSize', 20, 'LineWidth', 3)
hold on;
plot(real(Map(3,:)),imag(Map(3,:)), '-.','Color',[.3 .3 .3], 'MarkerSize', 20, 'LineWidth', 3)
end