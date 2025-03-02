classdef thin_airfoil_theory
    properties
    end
    methods
        function C_l = calculate_cl(obj, NACA, aoa, c)
            if nargin < 4
                c = 1;  % Default chord length
            end
            [A_0, A_1] = obj.calulate_A(NACA, aoa, c);

            C_l = 2*pi*(A_0 + A_1/2);
        end

        function [A_0, A_1] = calulate_A(obj, NACA, aoa, c)
            if nargin < 4
                c = 1;  % Default chord length
            end

            theta = [0, acos(1-2*NACA.p), pi];
        
            A_0_base = theta.*(2*NACA.p./c - 1) + sin(theta);
            A_0 = aoa ...
                - 1/pi*(NACA.m./NACA.p.^2.*A_0_base(2)-A_0_base(1) ...
                        + NACA.m./(1-NACA.p).^2.*A_0_base(3)-A_0_base(2));
            
            A_1 = obj.calulate_An(NACA, c, 1);
        end

        function A_n =  calulate_An(obj, NACA, c, n)
            theta = [0, acos(1-2*NACA.p), pi];
            base = sin((n+1).*theta)./(2.*(n+1)) ...
                    + 2*NACA.p.*sin(n.*theta)./(c.*n)...
                    - sin(n.*theta)./n...
                    + sin((1-n).*theta)./(2*(1-n));
            A_n = 1/pi*(NACA.m./NACA.p.^2.*(base(2)-base(1)) ...
                        + NACA.m./(1-NACA.p).^2.*(base(3)-base(2)));
        end

        function Cl = calculate_cl_numeric(obj, alpha, m, p, c)
            % alpha: angle of attack in radians
            % m: maximum camber (as a fraction of chord length)
            % p: location of maximum camber (as a fraction of chord length)
            % c: chord length
        
            % Number of points for numerical integration
            N = 1000;
            theta = linspace(0, pi, N);
            dtheta = theta(2) - theta(1);
        
            % Transformation of x
            x = c / 2 * (1 - cos(theta));
        
            % Calculate the camber line y_c(x)
            y_c = zeros(1, N);
            for i = 1:N
                if x(i) <= p * c
                    y_c(i) = (m / p^2) * (2 * p * (x(i) / c) - (x(i) / c)^2);
                else
                    y_c(i) = (m / (1 - p)^2) * (1 - 2 * p + 2 * p * (x(i) / c) - (x(i) / c)^2);
                end
            end
        
            % Calculate the derivative dy_c/dx
            dy_cdx = gradient(y_c, x);
        
            % Calculate A0 and A1 using numerical integration
            A0 = alpha - (1 / pi) * trapz(theta, dy_cdx);
            A1 = (2 / pi) * trapz(theta, dy_cdx .* cos(theta));
        
            % Calculate the lift coefficient Cl
            Cl = 2 * pi * (A0 + A1 / 2);
        end
    end
end