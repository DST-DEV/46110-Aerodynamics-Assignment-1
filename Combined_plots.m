clc;clear;
%% Load Data
xfoil_res = load('xfoil_exports\XFOIL_results.mat').airfoils;
thin_res = load('xfoil_exports\thinAirfoilTheory.mat').thinAirfoilTheory;
airfoil_names = ["2312", "2324", "4412", "4424"];

%% User input
savefigs = true;
plot_C_l = true;
plot_dC_p = true;
plot_C_p = true;
plot_polars = true;

exp_fld = 'plots';

%% Plot settings
cols = ["#0072BD", "#D95319", "#EDB120", "#77AC30"];  % Colors of the lines
markers = ["none", "+", "o", "diamond"];  % Markers for the four methods
ms = [8, 8, 4.5, 6.5];  % Marker size for the plots of the four methods
lw = [1, 1, 1, 1];  % Linewidth for the lines of the four methods
ax_col = [0.2, 0.2, 0.2];  % Color of accented axes
ax_lw = 1.5;  % Line width of accented axes
fs = 12;  % Plot font size
fig_count = 0;

%% Preparation
% Create export directory if it doesn't exist
if ~exist(exp_fld, 'dir')
    mkdir(exp_fld);
end

%% Plot C_l vs alpha
if plot_C_l
    for i = 1:length(airfoil_names)
        % Find index of airfoil in structs
        i_thin = find(strcmp({thin_res.name}, airfoil_names(i)));
        i_xfoil = find(strcmp({xfoil_res.name}, airfoil_names(i)));
        
        % Create plot
        figure(i+fig_count);
        cla; hold on; grid on;
        colororder(cols);
        ax = gca;
            
        % Plot C_l curves 
        % Plot thin airfoil theory
        plt_cl_thin = plot(thin_res(i_thin).aoa,thin_res(i_thin).cl, ...
                           LineWidth=lw(1), Marker=markers(1), MarkerSize=ms(1));
    
        % Plot panel method
        plt_cl_panel = plot(0,0, ...
                            LineWidth=lw(2), Marker=markers(2), MarkerSize=ms(2));
    
        % % Plot Xfoil results (Interpolate bc alpha step is way to fine)
        % alpha_plot = -5:.25:15;
        % C_l_free = interp1(xfoil_res(i_xfoil).C_ld_free.alpha, ...
        %                    xfoil_res(i_xfoil).C_ld_free.C_l, ...
        %                    alpha_plot, ...
        %                    'linear', 'extrap');
        % C_l_fixed = interp1(xfoil_res(i_xfoil).C_ld_fixed.alpha, ...
        %                     xfoil_res(i_xfoil).C_ld_fixed.C_l, ...
        %                     alpha_plot, ...
        %                     'linear', 'extrap');
        % plt_cl_free = plot(alpha_plot, C_l_free, ...
        %                     LineWidth=lw(3), Marker=markers(3), MarkerSize=ms(3));
        % plt_cl_fixed = plot(alpha_plot, C_l_fixed, ...
        %                     Linewidth=lw(4), Marker=markers(4), MarkerSize=ms(4));

        % Plot Xfoil results (without interpolation)
        plt_cl_free = plot(xfoil_res(i_xfoil).C_ld_free.alpha, ...
                           xfoil_res(i_xfoil).C_ld_free.C_l, ...
                           LineWidth=lw(3), Marker=markers(3), MarkerSize=ms(3));
         plt_cl_fixed = plot(xfoil_res(i_xfoil).C_ld_fixed .alpha, ...
                           xfoil_res(i_xfoil).C_ld_fixed .C_l, ...
                           LineWidth=lw(4), Marker=markers(4), MarkerSize=ms(4));
        
        % Highlight x=0 and y=0 grid lines
        x_ax = xline(0, Color=ax_col, LineWidth=ax_lw); % Thick vertical line at x=0
        y_ax = yline(0, Color=ax_col, LineWidth=ax_lw); % Thick horizontal line at y=0
    
        %Order the plots 
        ax.Children = [plt_cl_thin; plt_cl_panel; plt_cl_free; plt_cl_fixed; x_ax; y_ax];
        hold off; 
    
        % Plot labels
        set(gcf,'Color','White');
        set(ax,'FontSize',fs);
        legend([plt_cl_thin, plt_cl_panel, plt_cl_free, plt_cl_fixed], ...
            {'Thin airfoil theory', 'Panel method', 'XFOIL (free transition)', ...
            'XFOIL (fixed transition)'}, 'Location', 'northwest')
        xlabel('AoA [deg]', 'Interpreter', 'latex');
        ylabel('$C_l$', 'Interpreter', 'latex');
        set(ax, 'TickLabelInterpreter', 'latex');
        
        ylim('auto');
        xticks(-10:2:15);
    
        % Save figure
        if savefigs
            exp_name = fullfile(exp_fld, ...
                        sprintf('C_l_vs_alpha_%s.pdf', airfoil_names(i)));
            exportgraphics(gcf, exp_name, 'ContentType', 'vector', ...
                'BackgroundColor', 'none', 'Resolution', 300);
        end
    end
else
    disp('C_l vs alpha not plotted')
end

fig_count = 4;

%% Plot ΔC_p over x/c
if plot_dC_p
    for i = 1:length(airfoil_names)
        % Find index of airfoil in structs
        i_thin = find(strcmp({thin_res.name}, airfoil_names(i)));
        i_xfoil = find(strcmp({xfoil_res.name}, airfoil_names(i)));
        
        % Create plot
        figure(i+fig_count);
        cla; hold on; grid on;
        colororder(cols);
        ax = gca;
        
        % Plot ΔC_p curves 
        % Plot thin airfoil theory
        plt_dcp_thin = plot(thin_res(i_thin).xc,thin_res(i_thin).dCp, ...
                           LineWidth=lw(1), Marker=markers(1), MarkerSize=ms(1));
    
        % Plot panel method
        plt_dcp_panel = plot(0,0, ...
                            LineWidth=lw(2), Marker=markers(2), MarkerSize=ms(2));
    
        % Plot Xfoil results (with interpolation to reduce number of points)
        % xc_plot = linspace(0,1,30);
        % dC_p_free = interp1(xfoil_res(i_xfoil).C_p_free.xc, ...
        %                    xfoil_res(i_xfoil).C_p_free.dC_p, ...
        %                    xc_plot, ...
        %                    'linear', 'extrap');
        % dC_p_fixed = interp1(xfoil_res(i_xfoil).C_p_fixed.xc, ...
        %                     xfoil_res(i_xfoil).C_p_fixed.dC_p, ...
        %                     xc_plot, ...
        %                     'linear', 'extrap');
        % plt_dcp_free = plot(xc_plot, dC_p_free, ...
        %                     LineWidth=lw(3), Marker=markers(3), MarkerSize=ms(3));
        % plt_dcp_fixed = plot(xc_plot, dC_p_fixed, ...
        %                     Linewidth=lw(4), Marker=markers(4), MarkerSize=ms(4));
        
        % Plot Xfoil results (without interpolation)
        plt_dcp_free = plot(xfoil_res(i_xfoil).C_p_free.xc, ...
                            xfoil_res(i_xfoil).C_p_free.dC_p, ...
                            LineWidth=lw(3), Marker=markers(3), MarkerSize=ms(3));
        plt_dcp_fixed = plot(xfoil_res(i_xfoil).C_p_fixed.xc, ...
                             xfoil_res(i_xfoil).C_p_fixed.dC_p, ...
                             Linewidth=lw(4), Marker=markers(4), MarkerSize=ms(4));
        
        % Highlight x=0 and y=0 grid lines
        x_ax = xline(0, Color=ax_col, LineWidth=ax_lw); % Thick vertical line at x=0
        y_ax = yline(0, Color=ax_col, LineWidth=ax_lw); % Thick horizontal line at y=0
    
        %Order the plots 
        ax.Children = [plt_dcp_thin; plt_dcp_panel; plt_dcp_free; plt_dcp_fixed; x_ax; y_ax];
        hold off; 
    
        % Plot labels
        set(gcf,'Color','White');
        set(ax,'FontSize',fs);
        legend([plt_dcp_thin, plt_dcp_panel, plt_dcp_free, plt_dcp_fixed], ...
            {'Thin airfoil theory', 'Panel method', 'XFOIL (free transition)', ...
            'XFOIL (fixed transition)'}, 'Location', 'northeast')
        xlabel('$x/c$', 'Interpreter', 'latex');
        ylabel('$\Delta C_p$', 'Interpreter', 'latex');
        set(ax, 'TickLabelInterpreter', 'latex');
        
        ylim('auto');
        xticks(0:.1:1);
    
        % Save figure
        if savefigs
            exp_name = fullfile(exp_fld, ...
                        sprintf('dC_p_vs_xc_%s.pdf', airfoil_names(i)));
            exportgraphics(gcf, exp_name, 'ContentType', 'vector', ...
                'BackgroundColor', 'none', 'Resolution', 300);
        end
    end
else
    disp('dC_p vs x/c not plotted')
end
fig_count = 8;

%% Plot C_p over x/c

if plot_C_p
    for i = 1:length(airfoil_names)
        % Find index of airfoil in structs
        i_xfoil = find(strcmp({xfoil_res.name}, airfoil_names(i)));
        
        % Create plot
        figure(i+fig_count);
        cla; hold on; grid on;
        colororder(cols(2:end));
        ax = gca;
        
        % Plot C_p curves 
        % Plot panel method
        plt_cp_panel = plot(0,0, ...
                            LineWidth=lw(2), Marker=markers(2), MarkerSize=ms(2));
    
        % Plot Xfoil results (with interpolation to reduce number of points)
        % xc_plot = linspace(0,1,30);
        % y_xfoil_free = xfoil_res(i_xfoil).C_p_free.y;
        % C_p_free_upper = interp1(xfoil_res(i_xfoil).C_p_free.x(y_xfoil_free>=0), ...
        %                    xfoil_res(i_xfoil).C_p_free.C_p(y_xfoil_free>=0), ...
        %                    xc_plot, ...
        %                    'linear', 'extrap');
        % C_p_free_lower = interp1(xfoil_res(i_xfoil).C_p_free.x(y_xfoil_free<=0), ...
        %                    xfoil_res(i_xfoil).C_p_free.C_p(y_xfoil_free<=0), ...
        %                    flip(xc_plot), ...
        %                    'linear', 'extrap');
        % C_p_free = horzcat(C_p_free_upper, C_p_free_lower);
        % 
        % y_xfoil_fixed = xfoil_res(i_xfoil).C_p_fixed.y;
        % C_p_fixed_upper = interp1(xfoil_res(i_xfoil).C_p_fixed.x(y_xfoil_fixed>=0), ...
        %                    xfoil_res(i_xfoil).C_p_fixed.C_p(y_xfoil_fixed>=0), ...
        %                    xc_plot, ...
        %                    'linear', 'extrap');
        % C_p_fixed_lower = interp1(xfoil_res(i_xfoil).C_p_fixed.x(y_xfoil_fixed<=0), ...
        %                    xfoil_res(i_xfoil).C_p_fixed.C_p(y_xfoil_fixed<=0), ...
        %                    flip(xc_plot), ...
        %                    'linear', 'extrap');
        % C_p_fixed = horzcat(C_p_fixed_upper, C_p_fixed_lower);
        % 
        % xc_plot = horzcat(xc_plot, flip(xc_plot));
        % plt_cp_free = plot(xc_plot, C_p_free, ...
        %                     LineWidth=lw(3), Marker=markers(3), MarkerSize=ms(3));
        % plt_cp_fixed = plot(xc_plot, C_p_fixed, ...
        %                     Linewidth=lw(4), Marker=markers(4), MarkerSize=ms(4));
        
        % Plot Xfoil results (without interpolation)
        plt_cp_free = plot(xfoil_res(i_xfoil).C_p_free.x, xfoil_res(i_xfoil).C_p_free.C_p, ...
                            LineWidth=lw(3), Marker=markers(3), MarkerSize=ms(3));
        plt_cp_fixed = plot(xfoil_res(i_xfoil).C_p_fixed.x, xfoil_res(i_xfoil).C_p_fixed.C_p, ...
                            Linewidth=lw(4), Marker=markers(4), MarkerSize=ms(4));
        
        % Highlight x=0 and y=0 grid lines
        x_ax = xline(0, Color=ax_col, LineWidth=ax_lw); % Thick vertical line at x=0
        y_ax = yline(0, Color=ax_col, LineWidth=ax_lw); % Thick horizontal line at y=0
    
        %Order the plots 
        ax.Children = [plt_cp_panel; plt_cp_free; plt_cp_fixed; x_ax; y_ax];
        hold off; 
    
        % Plot labels
        set(gcf,'Color','White');
        set(ax,'FontSize',fs);
        legend([plt_cp_panel, plt_cp_free, plt_cp_fixed], ...
            {'Panel method', 'XFOIL (free transition)', ...
            'XFOIL (fixed transition)'}, 'Location', 'southeast')
        xlabel('$x/c$', 'Interpreter', 'latex');
        ylabel('$C_p$', 'Interpreter', 'latex');
        set(ax, 'TickLabelInterpreter', 'latex');
        
        ylim('auto');
        xticks(0:.1:1);
    
        % Save figure
        if savefigs
            exp_name = fullfile(exp_fld, ...
                        sprintf('C_p_vs_xc_%s.pdf', airfoil_names(i)));
            exportgraphics(gcf, exp_name, 'ContentType', 'vector', ...
                'BackgroundColor', 'none', 'Resolution', 300);
        end
    end
else
    disp('C_p vs x/c not plotted')
end
fig_count = 12;

%% Plot C_l vs C_d

if plot_polars
    for i = 1:length(airfoil_names)
        % Find index of airfoil in structs
        i_xfoil = find(strcmp({xfoil_res.name}, airfoil_names(i)));
        
        % Create plot
        figure(i+fig_count);
        cla; hold on; grid on;
        colororder(cols(3:end));
        ax = gca;
        
        % Plot Xfoil results
        plt_pol_free = plot(xfoil_res(i_xfoil).C_ld_free.C_d, xfoil_res(i_xfoil).C_ld_free.C_l, ...
                            LineWidth=lw(3), Marker=markers(3), MarkerSize=ms(3));
        plt_pol_fixed = plot(xfoil_res(i_xfoil).C_ld_fixed.C_d, xfoil_res(i_xfoil).C_ld_fixed.C_l, ...
                            Linewidth=lw(4), Marker=markers(4), MarkerSize=ms(4));
        
        % Highlight x=0 and y=0 grid lines
        x_ax = xline(0, Color=ax_col, LineWidth=ax_lw); % Thick vertical line at x=0
        y_ax = yline(0, Color=ax_col, LineWidth=ax_lw); % Thick horizontal line at y=0
    
        %Order the plots 
        ax.Children = [plt_pol_free; plt_pol_fixed; x_ax; y_ax];
        hold off; 
    
        % Plot labels
        set(gcf,'Color','White');
        set(ax,'FontSize',fs);
        legend([plt_pol_free, plt_pol_fixed], ...
            {'Free transition', 'Fixed transition'}, 'Location', 'southeast')
        xlabel('$C_d$', 'Interpreter', 'latex');
        ylabel('$C_l$', 'Interpreter', 'latex');
        set(ax, 'TickLabelInterpreter', 'latex');
        
        ylim('auto');
    
        % Save figure
        if savefigs
            exp_name = fullfile(exp_fld, ...
                        sprintf('C_l_vs_C_d_%s.pdf', airfoil_names(i)));
            exportgraphics(gcf, exp_name, 'ContentType', 'vector', ...
                'BackgroundColor', 'none', 'Resolution', 300);
        end
    end
else
    disp('C_l vs C_d not plotted')
end
fig_count = 16;