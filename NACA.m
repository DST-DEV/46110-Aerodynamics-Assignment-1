clear all; close all; clc;

function y = NACA(m, p, c)
    
    x = linspace(0,1,500)

    y_1 = m./p.^2 * (2*p.*(x./c) - (x./c).^2);
    y_2 = m./(1-p).^2 * (1 - 2*p + 2*p.*x./c - (x./c).^2);
    y = horzcat(y_1(x<=p), y_2(x>=p));
end

    