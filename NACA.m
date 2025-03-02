
classdef NACA
    properties
      m
      p
      t 
      N = 200
    end
    methods
        %% --- Class constructor
        function obj = NACA(m, p, t, N)
            obj.m = m;
            obj.p = p;
            obj.t = t;
            if nargin == 4
                obj.N = N;
            end
        end

        %% --- Function to calculate the outer shape and camber line of a 4-digit NACA ---
        function [x, yc, xu, yu, xl, yl] = naca_airfoil(obj, c, N)
            if nargin < 3
                N = obj.N;  % Default number of points
            end
            if nargin < 2
                c = 1;  % Default chord length
            end
            x = linspace(0, 1, N)*c;
        
            % Get camber line and local slope
            [yc, theta] = obj.naca_camber(x, c);
        
            % Get thickness distribution
            yt = obj.naca_thickness(x);
        
            % Upper surface
            xu = x - yt .* sin(theta);
            yu = yc + yt .* cos(theta);
        
            % Lower surface
            xl = x + yt .* sin(theta);
            yl = yc - yt .* cos(theta);
        end
            
        %% --- Function to the camber line and slope ---
        % naca_camber(x, m, p) returns:
        %   yc:   camber-line y-coordinate
        %   theta: local slope angle = atan(dy/dx)
        %
        % x, m, p are dimensionless (0 <= xc <= 1).
        function [yc, theta] = naca_camber(obj, x, c)
            % Preallocate
            yc    = zeros(size(x));
            dycDx = zeros(size(x));
            
            xc = x./c;

            % Indices for piecewise definition
            idx1 = (xc < obj.p);
            idx2 = (xc >= obj.p);
        
            % For 0 <= x < p
            yc(idx1)    = (obj.m ./ obj.p^2) ...
                .* (2*obj.p*xc(idx1) - xc(idx1).^2);
            dycDx(idx1) = (2*obj.m ./ obj.p^2) .* (obj.p - xc(idx1));
        
            % For p <= x <= 1
            yc(idx2)    = (obj.m ./ (1 - obj.p)^2) ...
                 .* (1 - 2*obj.p + 2*obj.p*xc(idx2) - xc(idx2).^2);
            dycDx(idx2) = (2*obj.m ./ (1 - obj.p)^2) .* (obj.p - xc(idx2));
        
            % Slope angle
            theta = atan(dycDx);
        end
        
        %% --- Function to calculate the thickness distribution function for 4-digit NACA ---
        % naca_thickness(x, t) returns:
        %   yt: thickness distribution at each x
        %
        % x, t are dimensionless.
        function yt = naca_thickness(obj, x)
            yt = 5 * obj.t .* ( ...
                0.2969 * sqrt(x) ...
              - 0.1260 * x ...
              - 0.3516 * x.^2 ...
              + 0.2843 * x.^3 ...
              - 0.1015 * x.^4 );
        end
    end
end