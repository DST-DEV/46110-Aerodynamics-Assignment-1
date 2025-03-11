classdef XFOIL_NACA
    properties
        NACA
        xfoil_exe
        Re = 1.5e6;
        Ma = 0;
        n_crit = 9;
        Xt = .1;
    end
    methods
        %% --- Class constructor
        function obj = XFOIL_NACA(NACA, xfoil_exe)
            obj.NACA = NACA;
            if nargin<2
                obj.xfoil_exe = fullfile('..','..', ...
                    '03_additional_material', 'Xfoil', 'xfoil.exe');  %Path to Xfoil executable
            else
                obj.xfoil_exe = xfoil_exe;
            end
        end
        
        %% --- Function to calculate C_l and C_d of 4-digit NACA for a range of angles of attack---
        function [alpha, C_l, C_d] = calc_cl (obj, free_BL, use_cached_input, ...
                AoAs, numNodes, max_iter)
            if nargin < 6
                max_iter = 200;
            end
            if nargin < 5
                numNodes = 500;
            end
            if nargin < 4
                AoAs = [-5, 15, .25];  % Angles of attack [deg]
            end
            if nargin < 3
                use_cached_input = false;
            end
            if nargin < 2
                free_BL = false;
            end
            
            % Xfoil preparation & calculation
            exp_fld = 'xfoil_exports' 
            if free_BL
                fname_input = ['input_NACA' obj.NACA '_free_BL.txt'];  % XFoil input filename
                fname_cl_res = ['Cl_NACA' obj.NACA '_free_BL.txt'];  % Lift coefficient filename
            else
                fname_input = ['input_NACA' obj.NACA '_fixed_BL.txt'];  % Airfoil coordinates filename
                fname_cl_res = ['Cl_NACA' obj.NACA '_fixed_BL.txt'];  % Lift coefficient filename
            end
            
            fpath_coords = ['Coords_NACA' obj.NACA '.txt'];  % Airfoil coordinates filename
            
            % Create export directory if it doesn't exist
            if ~exist(exp_fld, 'dir')
                mkdir(exp_fld);
            end

            % Delete files if they exist
            if (exist(fullfile(exp_fld, fpath_coords),'file'))
                delete(fullfile(exp_fld, fpath_coords));
            end
            if (exist(fullfile(exp_fld, fname_cl_res),'file'))
                delete(fullfile(exp_fld, fname_cl_res));
            end

            if use_cached_input
                assert(exist(fullfile(exp_fld, fname_input),'file'), "Input file not found")
            else
                % Prepare input to xfoil
                % Create the airfoil
                fid = fopen(fullfile(exp_fld, fname_input), 'w');
                fprintf(fid, ['NACA ' obj.NACA '\t\t\t\t\t! Airfoil selection\n']);
                fprintf(fid, 'PPAR\t\t\t\t\t\t! show paneling menu\n');
                fprintf(fid, sprintf('N %d\t\t\t\t\t\t! Number of panels\n', numNodes));
                fprintf(fid, '\n\n');
                
                % Save the airfoil data points
                fprintf(fid,['PSAV ' exp_fld '\\' fpath_coords '\n']);
                
                % Calculate C_l vs alpha values
                fprintf(fid, 'OPER\t\t\t\t\t\t! Enter operational menu\n');
                fprintf(fid, sprintf('ITER %d\t\t\t\t\t\t\t! Maximum no. of iterations\n', max_iter));
                fprintf(fid, sprintf('VISC %10e\t\t\t\t! Set Reynolds number\n', obj.Re));
                fprintf(fid, sprintf('Mach %.4f\t\t\t\t\t! Set Mach number\n', obj.Ma));
                
                fprintf(fid,'VPAR\t\t\t\t\t\t! Enter BL parameter menu\n');
                if free_BL
                    fprintf(fid, sprintf(['N %d\t\t\t\t\t\t\t! '...
                        'Set critical amplification exponent (free transition)\n'], obj.n_crit));
                else
                    fprintf(fid, 'XTR\n');
                    fprintf(fid, sprintf('%.3f\n', obj.Xt));  % Set upper trip positions
                    fprintf(fid, '\n');
                end
                fprintf(fid, '\n');
                
                % Configure output file
                fprintf(fid, 'PACC\t\t\t\t\t\t! Enable polar accumulation (stores results)\n');
                fprintf(fid, [exp_fld '\\' fname_cl_res '\n']);
                fprintf(fid, '\n');
                
                %Run calculation for all AoAs
                fprintf(fid, sprintf('ASEQ %s\n' , sprintf('%.2f ' , AoAs)));
                
                %Finish up
                fprintf(fid, 'PACC\t\t\t\t\t\t! Disable polar accumulation (closes output file)\n');
                fprintf(fid, '\nQUIT\n');
                
                % Close file
                fclose(fid);

                % Run XFoil using input file
                cmd = obj.xfoil_exe + " < " + fullfile(exp_fld, fname_input);
                %[status,result] = system(cmd);
                system(cmd)

                [alpha, C_l, C_d] = obj.read_C_ld (free_BL);
            end
        end

        %% --- Function to load the xy coordinates of a 4-digit NACA from a text file from Xfoil---
        function [XB, YB] = read_coords (obj)
            % Determine filename
            fpath_coords = ['xfoil_exports\Coords_NACA' obj.NACA '.txt'];

            % Read airfoil coordinates
            coords_airfoil_file = fopen(fpath_coords);

            dataBuffer = textscan(coords_airfoil_file,'%f %f','CollectOutput',1,...
                                             'Delimiter','','HeaderLines',0);
            fclose(coords_airfoil_file);

            % Separate boundary points
            XB = dataBuffer{1}(:,1);
            YB = dataBuffer{1}(:,2); 
        end

        %% --- Function to load the AoAs, C_l and C_d values of a 4-digit NACA from a text file from Xfoil---
        function [alpha, C_l, C_d] = read_C_ld (obj, free_BL)
            % Determine filename
            if free_BL
                fpath_res = ['xfoil_exports\Cl_NACA' obj.NACA '_free_BL.txt'];
            else
                fpath_res = ['xfoil_exports\Cl_NACA' obj.NACA '_fixed_BL.txt'];
            end

            % Read lift coefficients
            res_file = fopen(fpath_res);
            dataBuffer = textscan(res_file,'%f %f %f %f %f %f %f', ...
                                        'HeaderLines',12,...
                                        'CollectOutput',1,...
                                        'Delimiter','');
            fclose(res_file);
            
            % Separate Cp data
            alpha  = dataBuffer{1,1}(:,1); 
            C_l  = dataBuffer{1,1}(:,2); 
            C_d = dataBuffer{1,1}(:,3);   % Boundary point Y-coordinate
        end
        
        %% --- Function to plot the shape of an airfoil---
        function plot_airfoil (obj, XB, YB)
            % Split airfoil into (U)pper and (L)ower
            XB_U = XB(YB >= 0);
            XB_L = XB(YB < 0);
            YB_U = YB(YB >= 0);
            YB_L = YB(YB < 0);
            
            % Plot: Airfoil
            figure(1);
            cla; hold on; grid off;
            set(gcf,'Color','White');
            set(gca,'FontSize',12);
            plot(XB_U,YB_U,'b.-');
            plot(XB_L,YB_L,'r.-');
            xlabel('X Coordinate');
            ylabel('Y Coordinate');
            axis equal;
        end
        
        %% --- Function to plot the lift coefficient of an airfoil over a range of AoAs---
        function plot_Cl (obj, alpha, C_l)
            % Plot lift coefficient
            figure(2);
            cla; hold on; grid on;
            set(gcf,'Color','White');
            set(gca,'FontSize',12);
            
            plot(alpha , C_l, 'bo-', 'LineWidth', 2, 'MarkerSize', 6);
            
            % Highlight x=0 and y=0 grid lines
            gray_color = [0.2, 0.2, 0.2];
            xline(0, 'Color', gray_color, 'LineWidth', 1.5); % Thick vertical line at x=0
            yline(0, 'Color', gray_color, 'LineWidth', 1.5); % Thick horizontal line at y=0
            
            % Plot labels
            xlabel('AoA [deg]', 'Interpreter', 'latex');
            ylabel('$C_l$', 'Interpreter', 'latex');
            set(gca, 'TickLabelInterpreter', 'latex');
            
            ylim('auto');
            xticks(-10:2:15);
        end
    end
end