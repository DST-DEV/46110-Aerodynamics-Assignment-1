%% ASSIGMENT 1 46110: AIRFOIL CHARACTERISTICS USING DIFFERENT ANALYSIS METHODS

clear
clc
close all


%% DATA: NACA 2312 // NACA 2324 // NACA 4412 // NACA 4424

num = [1 2 3 4]; % Value to change the airfoil data
n_div = 50; % Number of divisions for the arrays
check_aoa = 10; % Angle of attack at which we want to check the deltaCp

c = 1; % [-] Chord line of the airfoils
c_t = pi; % [rad] Transformed chord line
density = 1.225; % [kg/m3] Air density
Qinf = 1; % [m/s] Free stream velocity

x = linspace(0,1,n_div); % Array to calculate the camber line
theta = acos(1 - 2.*x./c);
aoa_deg = linspace(-10,15,n_div);
aoa = deg2rad(aoa_deg);

%% THIN AIRFOIL THEORY

for i = 1:length(num) % We loop to calculate the values for the four different airfoils

    % Function to get the airfoil data
    [m,p,t] = airfoil_selection(num(i));
    
    % Calculating the transformed derivative of the camber line
    dycdx_t = camber_transformed_derivation(m,p,c_t,theta);
     
    % Calculating the fourier components for Glauert's solution
    [A0, A1] = fourier_components(dycdx_t,aoa,theta);
    
    % Calculating the lift coefficient of the airfoil by Glauert's solution
    cl(i,:) = lift_coefficient(A0, A1);
    
    % Finding the index of the angle of attack equivalent to 10 degrees
    [~,ind_aoa] = min(abs(rad2deg(aoa)-check_aoa));
    
    % Calculating the variation in the pressure coefficient
    deltaCp(i,:) = delta_pressure_coefficient(dycdx_t, aoa(ind_aoa), Qinf, theta);
end
    
%% SAVING THE RESULTS

% Creating a struct to save the data of all four airfoils as one variable

thinAirfoilTheory = [struct('name', '2312', 'xc', x, 'aoa', aoa_deg, 'cl', cl(1,:), 'dCp',deltaCp(1,:));
                     struct('name', '2324', 'xc', x, 'aoa', aoa_deg, 'cl', cl(2,:), 'dCp',deltaCp(2,:));
                     struct('name', '4412', 'xc', x, 'aoa', aoa_deg, 'cl', cl(3,:), 'dCp',deltaCp(3,:));
                     struct('name', '4424', 'xc', x, 'aoa', aoa_deg, 'cl', cl(4,:), 'dCp',deltaCp(4,:))];


%% PLOTTING THE RESULTS

figure(1)
plot(rad2deg(aoa),cl,'k',LineWidth=2)
grid on
xlabel('Angle of attack [deg]')
ylabel('Lift coefficient [-]')

figure(2)
plot(x,deltaCp,'k',LineWidth=2)
grid on
xlabel('x/c')
ylabel('\delta Cp')











%% FUNCTIONS

function [m,p,t] = airfoil_selection(number)
    %{
    Returns the airfoil data to analyse depending on the choice of the
    user; 1 = NACA 2312, 2 = NACA 2324, 3 = NACA 4412 and 4 = NACA 4424.

    Args:
        number (str): Value between 1 and 4 that equals each airfoil type

    Raises:
        An error if the value is not an integer
        An error if the value is not between 1 and 4

    Returns:
        m (float): The maximum camber of the airfoil, as a percentage of
        chord
        p (float): Location of the maximum camber as a percentage of chord
        t (float): Maximum thickness of the airfoil as a percentage of
        chord
    %}

    % Array of the data from the four different airfoils
    max_camber = [0.02 0.02 0.04 0.04];
    location_max_camber = [0.3 0.3 0.40 0.4];
    max_thickness = [0.12 0.24 0.12 0.24];

    %Checking if the number is an integer and between 1 and 4
    if mod(number, 1) == 0
        if number > 0 && number <= 4
            m = max_camber(number);
            p = location_max_camber(number);
            t = max_thickness(number);
        else
            disp('The number must be between 1 and 4')
        end
    else
        disp('The number must be an integer')
    end

end

% ---------------------------------------------------------

% Function to calculate the transformed derivative of the camber line
function dycdx_t = camber_transformed_derivation(m, p, c_t, theta)
    %{
    Returns an array of values of the camber line equation derivated over
    theta, transformed into polar coordinates instead of cartesian as in:
    x = (1-cos(theta))*c/2.

    Args:
        m (float): The maximum camber of the airfoil, as a percentage of
        chord
        p (float): Location of the maximum camber as a percentage of chord
        c_t (float): Transformed value of the chord, equal to pi if x = 1
        theta (float): Transformed x coordinate

    Raises:
        Nothing is checked in this function

    Returns:
        dycdx_t (float): An array containing the theta derivative of the
        camber line equation
    %}

for jj = 2:length(theta)
        % Calculating the derivative based on the analitical expression
        if theta(jj) <= acos(abs(2*p/c_t - 1))
            dycdx_t(jj) = m/(p^2)*(2*p/c_t - (1-cos(theta(jj)))/c_t);
        else
            dycdx_t(jj) = m/((1-p)^2)*(2*p/c_t - (1-cos(theta(jj)))/c_t);
        end
end

end

% ---------------------------------------------------------

function [A0, A1] = fourier_components(dycdx_t, aoa, theta)
    %{
    Returns the Fourier components that compose Glauert's solution to the
    thin airfoil problem

    Args:
        dycdx_t (float): An array containing the theta derivative of the
        camber line equation
        aoa (float): The angle of attack at which we want to calculate the
        lift coefficient, in radians.
        theta (float): Transformed x coordinate

    Raises:
        Nothing is checked in this function

    Returns:
        A0 (float): First fourier component, dependent on the angle of
        attack
        A1 (float): Second fourier component
    %}
    
    % Calculating the Fourier components via integration
    A0 = aoa - (1/pi)*trapz(dycdx_t,theta);
    A1 = (2/pi)*trapz(dycdx_t.*cos(theta),theta);

end

% ---------------------------------------------------------

function cl = lift_coefficient(A0, A1)
    %{
    Returns the lift coefficient of the airfoil, calculated based on
    Glauert's solution and Fourier components

    Args:
        A0 (float): First fourier component, dependent on the angle of
        attack
        A1 (float): Second fourier component

    Raises:
        Nothing is checked in this function

    Returns:
        cl (float): Lift coefficient of the selected airfoil
    %}
    
    % Calculating the lift coefficient
    cl = 2*pi*(A0 + A1/2);
end

% ---------------------------------------------------------

function deltaCp = delta_pressure_coefficient(dycdx_t, aoa, Qinf, theta)
    %{
    Returns the variation in the pressure coefficient between the upper and
    lower part of the airfoil along the x axis in polar coordinates,  for a
    defined angle of attack.

    Args:
        dycdx_t (float): An array containing the theta derivative of the
        camber line equation
        ind_aoa (float): The angle of attack at which we want to calculate the
        lift coefficient, in radians.
        Qinf (float): Free stream wind speed
        theta (float): Transformed x coordinate       

    Raises:
        Nothing is checked in this function

    Returns:
        deltaCp (float): Array of the variation in the pressure coefficient
        between the upper and lower part of the airfoil
    %}

    % Calculating the Fourier components for the angle of attack selected
    [A0, A1] = fourier_components(dycdx_t,aoa,theta);

    % Calculating the Glauert's solution for the circulation distribution
    gamma = 2*Qinf*(A0.*((1 + cos(theta))./sin(theta)) + A1*sin(theta));

    % Calculating the variation in the pressure coefficient
    deltaCp = 2*gamma./Qinf;

end