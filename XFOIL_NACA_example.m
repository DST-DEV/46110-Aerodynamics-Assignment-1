clear; clc;
xfoil_exe = fullfile('..','..', ...
                    '03_additional_material', 'Xfoil', 'xfoil.exe');  %Path to Xfoil executable

airfoils = [ ...
    struct('name','2312', 'xfoil', XFOIL_NACA('2312', xfoil_exe)),
    struct('name','2324', 'xfoil', XFOIL_NACA('2324', xfoil_exe)),
    struct('name','4412', 'xfoil', XFOIL_NACA('4412', xfoil_exe)),
    struct('name','4424', 'xfoil', XFOIL_NACA('4424', xfoil_exe))
    ];

%% Run Calculations
rerun_calc = true;

AoAs = [-5, 15, .1];
numNodes = 500;
max_iter = 400;
use_cached_input = false;


if rerun_calc
    for i = 1:size(airfoils)
        % Run calculation with free transition
        [alpha, C_l, C_d] = airfoils(i).xfoil.calc_cl(true, use_cached_input, AoAs, numNodes, max_iter);
        airfoils(i).free_trans = struct('alpha', alpha, 'C_l', C_l, 'C_d', C_d);

        % Run calculation with fixed transition
        [alpha, C_l, C_d] = airfoils(i).xfoil.calc_cl(false, use_cached_input, AoAs, numNodes, max_iter);
        airfoils(i).fixed_trans = struct('alpha', alpha, 'C_l', C_l, 'C_d', C_d);
    end
end

