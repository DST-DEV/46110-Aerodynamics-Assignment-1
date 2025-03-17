clear; clc;

airfoils = [ ...
    struct('name','2312', 'NACA', NACA(0.02, 0.3, 0.12)),
    struct('name','2324', 'NACA', NACA(0.02, 0.3, 0.24)),
    struct('name','4412', 'NACA', NACA(0.04, 0.4, 0.12)),
    struct('name','4424', 'NACA', NACA(0.04, 0.4, 0.24))
    ];

%% Run Calculations
AoAs = -10:.1:15;
N_panels = 150;

% If output file already exists, load it
fpath_res = 'xfoil_exports\panel_results_davis.mat';

%% Calculate C_l & C_d
for i = 1:size(airfoils)
    % Run calculation with free transition
    [alpha, C_l, C_p, x, y] = calc_panel(airfoils(i).NACA, AoAs, N_panels);
    [xc, dC_p] = calc_dc_p (x, y, C_p(201,:));

    airfoils(i).panel_res = struct('alpha', alpha, 'C_l', C_l, ...
        'C_p', C_p, 'x', x, 'y', y, 'dC_p', dC_p, 'xc', xc);

end

save(fpath_res, 'airfoils');  % Save results

%% --- Function to calculate C_l and C_d of 4-digit NACA for a range of angles of attack---
function [alpha, C_l, C_p, xp, yp] = calc_panel (NACA, alpha, N_panels)
    % Get airfoil coordinates
     [x, yc, xu, yu, xl, yl, dycDx] = NACA.naca_airfoil(1, N_panels);
    x = [flip(xu), xl(2:end)];
    y = [flip(yu), yl(2:end)];
    n = length(x)-1;
    % Close the airfoil shape
    y (1) = 0;
    y (end) = 0;
    x(end) = 1;
    x(1) = 1;


    
    % Definition of panels and tangent vectors
    for j=1:n
        % Panel lengths
        plength(j) = sqrt((x(j+1)-x(j)).^2+(y(j+1)-y(j)).^2);
        
        % Control point (in the middle)
        xp(j)=0.5*(x(j+1)+x(j));
        yp(j)=0.5*(y(j+1)+y(j));
        
        % Tangent vectors
        Tx(j)=-(x(j+1)-x(j))./plength(j); 
        Ty(j)=-(y(j+1)-y(j))./plength(j);
        
        % Normal vectors
        Nx(j)=-Ty(j);
        Ny(j)= Tx(j);
    end
    
    % A and B matrices
    for i=1:n
        for j=1:n
           if i == j
               Ux(i,i)=0.0; % Self-induction eq. (19) in note
               Uy(i,i)=0.5; % Self-induction eq. (19) in note
           else
               % Induction at control point 'i' in local system 'j'
               % Eqs. (16) and (18) in note
               sx(i,j)=(xp(i)-xp(j))*Tx(j)+(yp(i)-yp(j))*Ty(j);
               sy(i,j)=(xp(i)-xp(j))*Nx(j)+(yp(i)-yp(j))*Ny(j);
               Ux1(i,j)=log( ((sx(i,j)+0.5*plength(j)).^2 + sy(i,j).^2)/((sx(i,j)-0.5*plength(j)).^2 + sy(i,j).^2) )/(4.*pi);
               Uy1(i,j)=( atan( (sx(i,j)+0.5*plength(j))/sy(i,j) ) - atan( (sx(i,j)-0.5*plength(j))/sy(i,j) ) )/(2.*pi);
              % Induction in global system
               Ux2(i,j)=Ux1(i,j)*Tx(j)-Uy1(i,j)*Ty(j);
               Uy2(i,j)=Ux1(i,j)*Ty(j)+Uy1(i,j)*Tx(j);
              % Induction in local system 'i'
               Ux(i,j)=Ux2(i,j)*Tx(i)+Uy2(i,j)*Ty(i);
               Uy(i,j)=Ux2(i,j)*Nx(i)+Uy2(i,j)*Ny(i);
           end
           % Eq. (24) in note by setting sigma=1
           A(i,j)=Uy(i,j); % Influence coefficients in A matrix (normal velocity)
           B(i,j)=Ux(i,j); % Influence coefficients in B matrix (tangential velocity)
        end
    end
    
    for i_a=1:length(alpha)
        alpha_i = deg2rad(alpha(i_a));
        % Boundary conditions
        F = -(Nx.*cos(alpha_i)+Ny.*sin(alpha_i));
        
        % Solution of system (solution of eq. (26))
        M = A\F';
        
        for i=1:n
            sum1=0.0;
            sum2=0.0;
            for j=1:n
                sum1=sum1+B(i,j)*M(j,1);
                sum2=sum2+A(i,j)*M(j,1);
            end
            Vt1(i)=sum1;
            Vt(i)=sum1+Tx(i)*cos(alpha_i)+Ty(i)*sin(alpha_i); %Tangential velocity (eq. 27)
            Vn(i)=sum2+Nx(i)*cos(alpha_i)+Ny(i)*sin(alpha_i); %Check of normal velocity (should be equal to zero)
        end
        
        %Vortex distribution
        sum=0.0;
        for i=1:n
            sum=sum+(i-1)*(n-i)*plength(i);
        end
        for i=1:n
            vort(i)=(i-1)*(n-i)/sum; %parabolic vortex distribution  (eq. 37)
        end
        
        C=A\B;
        D=B*C;
        Vrt=(A+D)*vort'; %Rotating onset flow (see eq. 40)
        Gamma = -(Vt1(1) + Vt1(end) ...
                + dot([cos(alpha_i), sin(alpha_i)], ...
                        ([Tx(1), Ty(1)] + [Tx(end), Ty(end)])))...
                ./ (Vrt(1) + Vrt(end));
        C_l(i_a) = 2.*Gamma/(abs(max(x)-min(x))); %Lift coefficient (K-J theorem: Cl=2*Gamma/Vc; V=1)
        Urt=Gamma*Vrt'+Vt; %Tangential velocity including circulation (eq. 41)
        for i=1:n
             C_p(i_a, i)=1.-Urt(i).^2; %Pressure coefficient with Kutta condition
        end
        
        % Verify Kutta condition
        assert(ismembertol(Urt(1), -Urt(end)), ...
            "Tangential velocities at trailing edge are not equal")
        assert (ismembertol(C_p(i_a, 1), C_p(i_a, end)), ...
            "Pressure coefficients don't match at trailing edge")
    end
end

%% --- Function to calculate ΔC_p of 4-digit NACA---
function [xc, dC_p] = calc_dc_p (x, y, C_p)
    % Split upper and lower surfaces
    upper_idx = y >= 0;
    lower_idx = y < 0;

    % Sort upper and lower surfaces by x-coordinates
    [x_upper, sort_idx_upper] = sort(x(upper_idx));
    C_p_upper = C_p(upper_idx);
    C_p_upper = C_p_upper(sort_idx_upper);
    
    [x_lower, sort_idx_lower] = sort(x(lower_idx));
    C_p_lower = C_p(lower_idx);
    C_p_lower = C_p_lower(sort_idx_lower);
    
    % Ensure x-coordinates match for interpolation
    C_p_lower_interp = interp1(x_lower, C_p_lower, x_upper, ...
        'linear', 'extrap');
    
    % Compute relative pressure coefficient (ΔCp)
    dC_p = C_p_lower_interp - C_p_upper;
    xc = x_upper;
end