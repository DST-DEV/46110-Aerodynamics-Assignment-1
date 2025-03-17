clear; clc; close all;
%% --- Main script portion: define airfoils and plot ---
% (m, p, t) all in fraction-of-chord terms

airfoilDefs = { ...
    'NACA 2312', NACA(0.02, 0.3, 0.12); ...
    'NACA 2324', NACA(0.02, 0.3, 0.24); ...
    'NACA 4412', NACA(0.04, 0.4, 0.12); ...
    'NACA 4424', NACA(0.04, 0.4, 0.24); ...
    };

colors = ["#0072BD", "#D95319", "#EDB120", "#77AC30"];  % Colors of the lines
fs = 16;  % Plot font size

figure('Name','NACA 4-Digit Airfoils','NumberTitle','off');
hold on; grid on;

for i = 1:size(airfoilDefs,1)
    name = airfoilDefs{i,1};

    [x, yc, xu, yu, xl, yl, dycDx] = airfoilDefs{i,2}.naca_airfoil();

    % Plot upper surface
    plot(xu, yu, 'Color', colors(i), 'LineWidth', 1.2, ...
        'DisplayName',[name]);
    % Plot lower surface
    plot(xl, yl, 'Color', colors(i), 'LineWidth', 1.2, ...
        'LineStyle','-', 'HandleVisibility','off');
    % Plot camber line (dashed)
    plot(x, yc, '--', 'Color', colors(i), ...
        'LineWidth', 1.0, 'HandleVisibility','off');
end

% Configure ticks
xticks(0:.1:1)
yticks(-.3:.1:.3)

% Configure labels and legend
set(gca,'FontSize',fs);
xlabel('$x/c$', 'Interpreter', 'latex');
ylabel('$y/c$', 'Interpreter', 'latex');
axis equal;
legend('Location','southeast', 'Interpreter', 'latex', ...
       'NumColumns', 2);
set(gca, 'TickLabelInterpreter', 'latex');
hold off;
exportgraphics(gcf, 'plots/NACA_comparison.pdf', 'ContentType', 'vector', ...
                'BackgroundColor', 'none', 'Resolution', 300);