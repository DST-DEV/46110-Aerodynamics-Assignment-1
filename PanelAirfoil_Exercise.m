% Panel method for Exersize 5
% Jens N. Sørensen: Updated 27-02-2025
% The equations are based on the note 'Panel Methods' by Jens N. Sørensen
clear all

% Angle of attack
alpha_deg = 0.; % Angle of attack in degrees
alpha = alpha_deg*pi/180;
Ix=1; %Ix=0: Compute circle; Ix=1: Read data from airfoil
if Ix == 1
    % Input of airfoil contour
    %xdat = load('KT_airfoil.txt');[n1,n2]=size(xdat);
    xdat = load('NACA0018.txt');
    [n1,n2]=size(xdat);
    x=xdat(:,1);
    y=xdat(:,2);
    n=n1-1;
else
    % Generation of circle contour; t goes from 0 to 2*pi
    n1=300;  % Number of points
    n=n1-1; % Number of panels
    dt=2*pi/n;
    t(1)=0.;
    x(1)=1.;
    y(1)=0.;
    for i=2:n1
        t(i)=dt*(i-1);
        x(i)=cos(t(i));
        y(i)=sin(t(i));
    end
end

% Plot of airfoil (alternative, a circle)
ifig=1;
figure(ifig)
plot(x,y,'-ko','LineWidth',2,'MarkerSize',6);
axis equal
title('Geometry');
grid
hold off

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

% Boundary conditions
F = -(Nx.*cos(alpha)+Ny.*sin(alpha));

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
    Vt(i)=sum1+Tx(i)*cos(alpha)+Ty(i)*sin(alpha); %Tangential velocity (eq. 27)
    Vn(i)=sum2+Nx(i)*cos(alpha)+Ny(i)*sin(alpha); %Check of normal velocity (should be equal to zero)
    Cp(i)=1.-Vt(i).^2; %Pressure coefficient without Kutta condition 
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
        + dot([cos(alpha), sin(alpha)], ...
                ([Tx(1), Ty(1)] + [Tx(end), Ty(end)])))...
        ./ (Vrt(1) + Vrt(end));
Cl=2.*Gamma/(abs(max(x)-min(x))); %Lift coefficient (K-J theorem: Cl=2*Gamma/Vc; V=1)
Urt=Gamma*Vrt'+Vt; %Tangential velocity including circulation (eq. 41)
for i=1:n
     Cpr(i)=1.-Urt(i).^2; %Pressure coefficient with Kutta condition
end

ifig=ifig+1;
figure(ifig)
plot(xp,Vt,'-bo','LineWidth',2,'MarkerSize',6);
title('Velocity distribution without circulation');
grid
hold off

ifig=ifig+1;
figure(ifig)
plot(xp,-Cp,'-bo','LineWidth',2,'MarkerSize',6);
title('Pressure coefficient without circulation');
xlabel('x'); ylabel('-Cp');
grid
hold off

ifig=ifig+1;
figure(ifig)
plot(xp,Urt,'-bo','LineWidth',2,'MarkerSize',6);
title('Velocity distribution with circulation');
grid
hold off

ifig=ifig+1;
figure(ifig)
plot(xp,-Cpr,'-bo','LineWidth',2,'MarkerSize',6);
title('Pressure coefficient with circulation');
axis([min(x) max(x) -1 max(-Cpr)]);
xlabel('x'); ylabel('-Cp');
grid
hold off
