%% ASSIGNMENT 1 46110: AIRFOIL CHARACTERISTICS USING DIFFERENT ANALYSIS METHODS

clear
clc
close all

%% DATA: NACA 2312 // NACA 2324 // NACA 4412 // NACA 4424

num = [1 2 3 4]; % Airfoil selection numbers
n = 50; % Number of panels
c = 1; % Chord length [m]
aoa_deg = linspace(-10,15,n); % Angle of attack range [deg]
aoa = deg2rad(aoa_deg); % Convert to radians

%% PANEL METHOD

for i = 1:length(num) 
    [m, p, t] = airfoil_selection(num(i)); % Get airfoil parameters
    [x, y] = naca_airfoil(m, p, t, c, n); % Generate airfoil geometry
    cl(i,:) = panel_lift_coefficient(x, y, aoa); % Compute Cl using panel method
end

%% PLOTTING RESULTS

figure(1)
plot(rad2deg(aoa), cl, 'k', LineWidth=2)
grid on
xlabel('Angle of Attack [deg]')
ylabel('Lift Coefficient [-]')
legend('NACA 2312', 'NACA 2324', 'NACA 4412', 'NACA 4424', Location='best')

%% FUNCTIONS

function [m, p, t] = airfoil_selection(number)
    % Returns the airfoil parameters (m, p, t) for given NACA 4-digit code
    max_camber = [0.02 0.02 0.04 0.04];
    location_max_camber = [0.3 0.3 0.4 0.4];
    max_thickness = [0.12 0.24 0.12 0.24];

    if number >= 1 && number <= 4 && mod(number,1) == 0
        m = max_camber(number);
        p = location_max_camber(number);
        t = max_thickness(number);
    else
        error('Invalid airfoil selection. Choose a number between 1 and 4.')
    end
end

function [x, y] = naca_airfoil(m, p, t, c, n)
    % Generates NACA 4-digit airfoil coordinates
    
    x = linspace(0, c, n);
    yt = (t/0.2) * (0.2969*sqrt(x) - 0.1260*x - 0.3516*x.^2 + 0.2843*x.^3 - 0.1015*x.^4);

    yc = zeros(size(x));
    dyc_dx = zeros(size(x));

    for i = 1:length(x)
        if x(i) < p*c
            yc(i) = (m/p^2) * (2*p*x(i)/c - (x(i)/c)^2);
            dyc_dx(i) = (2*m/p^2) * (p - x(i)/c);
        else
            yc(i) = (m/(1-p)^2) * ((1-2*p) + 2*p*x(i)/c - (x(i)/c)^2);
            dyc_dx(i) = (2*m/(1-p)^2) * (p - x(i)/c);
        end
    end

    theta = atan(dyc_dx);
    xu = x - yt .* sin(theta);
    xl = x + yt .* sin(theta);
    yu = yc + yt .* cos(theta);
    yl = yc - yt .* cos(theta);

    x = [flip(xu), xl(2:end)];
    y = [flip(yu), yl(2:end)];
end

function cl = panel_lift_coefficient(x, y, aoa)
    % Computes lift coefficient using the panel method
    
    n = length(x) - 1;

    for j = 1:n
        plength(j) = sqrt((x(j+1) - x(j))^2 + (y(j+1) - y(j))^2);
        xp(j) = 0.5 * (x(j+1) + x(j));
        yp(j) = 0.5 * (y(j+1) + y(j));
        Tx(j) = -(x(j+1) - x(j)) / plength(j);
        Ty(j) = -(y(j+1) - y(j)) / plength(j);
        Nx(j) = -Ty(j);
        Ny(j) = Tx(j);
    end

    A = zeros(n, n);
    B = zeros(n, n);

    for i = 1:n
        for j = 1:n
            if i == j
                A(i, j) = 0.5;
                B(i, j) = 0;
            else
                sx = (xp(i) - xp(j)) * Tx(j) + (yp(i) - yp(j)) * Ty(j);
                sy = (xp(i) - xp(j)) * Nx(j) + (yp(i) - yp(j)) * Ny(j);
                Ux1 = log(((sx + 0.5 * plength(j))^2 + sy^2) / ((sx - 0.5 * plength(j))^2 + sy^2)) / (4 * pi);
                Uy1 = (atan((sx + 0.5 * plength(j)) / sy) - atan((sx - 0.5 * plength(j)) / sy)) / (2 * pi);
                Ux2 = Ux1 * Tx(j) - Uy1 * Ty(j);
                Uy2 = Ux1 * Ty(j) + Uy1 * Tx(j);
                Ux(i, j) = Ux2 * Tx(i) + Uy2 * Ty(i);
                Uy(i, j) = Ux2 * Nx(i) + Uy2 * Ny(i);
                A(i, j) = Uy(i, j);
                B(i, j) = Ux(i, j);
            end
        end
    end

    cl = zeros(1, length(aoa));

    for a = 1:length(aoa)
        alpha = aoa(a);
        F = -(Nx .* cos(alpha) + Ny .* sin(alpha));
        M = A \ F';

        for i = 1:n
            Vt(i) = sum(B(i, :) .* M') + Tx(i) * cos(alpha) + Ty(i) * sin(alpha);
        end

        sumVort = sum((0:n-1) .* (n-1:-1:0) .* plength);
        vort = ((0:n-1) .* (n-1:-1:0)) / sumVort;

        C = A \ B;
        D = B * C;
        Vrt = (A + D) * vort';
        Gamma = -(Vt(1) + Vt(end) + dot([cos(alpha), sin(alpha)], ([Tx(1), Ty(1)] + [Tx(end), Ty(end)]))) / (Vrt(1) + Vrt(end));

        cl(a) = 2 * Gamma / (max(x) - min(x));
    end
end
