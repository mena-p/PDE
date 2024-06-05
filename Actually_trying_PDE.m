clear variables;
close all;

% Intervalle
Intx = [200 500];
Inty = [-100 100];

% Dirichlet und Neumann Funktionen
gD = @(x,y) y*(1-y);
gN = @(x,y) -400;
Boundary_condition = @(x,y) 400;

% Quellterm
f = @(x,y) 0; 

n = 7;
nx = 2^n;
ny = 2^n;

% Gitterweite
hx = (Intx(2)-Intx(1))/nx;
hy = (Inty(2)-Inty(1))/ny;

% Aufstellen des Knotenvektors (einschließlich der Randpunkte)
Nx = nx+1;
Ny = ny+1;
nodes_x = linspace(Intx(1),Intx(2),Nx);
nodes_y = linspace(Inty(1),Inty(2),Ny);

% Diffusivity of water vapor in air
D = 0.000025;

% Updraft velocity vector (only vertical component x)
vx = zeros(Nx*Ny,1);
for i = 1:Nx
    for j = 1:Ny
        % coordinate where velocity is calculated
        x = nodes_x(i);
        y = nodes_y(j);

        % lexicographic index
        l = (j-1)*Nx+i;
        vx(l) = 1*allen(y,x,4,1000,250000,0); % y is location, x is height in our PDE!
    end
end
% arrumar a ordem lexicográfica movendo o array todo um para a esquerda (e eliminando o primeiro elemento no processo)
vx(1:Nx*Ny-1) = vx(2:Nx*Ny);

% remover NaN
vx(isnan(vx)) = 0;

%%% 5-Punkte Stern %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Aufstellen der Steifigkeitsmatrix
Np = Nx*Ny;
e = ones(Np,1);

A =  spdiags([D*hx*hx*e D*hy*hy*e+hy*hy*hx.*vx -2*D*(hy*hy+hx*hx)*e-hx*hy*hy*vx D*hy*hy*e D*hx*hx*e],...
            [-Nx,-1,0,1,Nx],Np,Np);
%%% Advektiver Teil: Zentrierte Finite-Differenzen
% a = 1;
% b = 12;
% B = spdiags([-0.5*b*hy*hx*hx*e -0.5*a*hx*hy*hy*e 0.5*a*hx*hy*hy*e 0.5*b*hy*hx*hx*e],...
%             [-Nx,-1,1,Nx],Np,Np);
% A = A + B;

% Aufstellen der rechten Seite
F = zeros(Np,1);
for j = 1:Ny
    for i = 1:Nx
        F((j-1)*Nx+i) = hx*hx*hy*hy*f((i-1)*hx,(j-1)*hy); 
    end
end

% Einsetzen der Neumann-Randbedingungen
for j=2:Ny-1 % Right edge of domain (top of thermal)
    A((j-1)*Nx+Nx,:) = 0;
    % Backwards difference to approximate derivative
    A((j-1)*Nx+Nx,(j-1)*Nx+Nx-1) = -1;
    A((j-1)*Nx+Nx,(j-1)*Nx+Nx) = 1;
    F((j-1)*Nx+Nx) = 0;%-hx*BC_discontinuous_FD((Nx-1)*hx,(j-1)*hy);
end
for j=2:Ny-1 % Left edge of domain (surface)
    A((j-1)*Nx+1,:) = 0;
    % Forwards difference to approximate derivative
    A((j-1)*Nx+1,(j-1)*Nx+1) = -1;
    A((j-1)*Nx+1,(j-1)*Nx+1+1) = 1;
    F((j-1)*Nx+1) = -hx*BC_discontinuous_FD((1-1)*hx,(j-1)*hy);
end

% Einsetzen der Dirichlet-Randbedingungen
for i = 1:Nx % Untere Kante
    A(i,:) = 0;
    A(i,i) = 1;
    F(i) = 0;
end
for i = 1:Nx % Obere Kante
    A((Ny-1)*Nx+i,:) = 0;
    A((Ny-1)*Nx+i,(Ny-1)*Nx+i) = 1;
    F((Ny-1)*Nx+i) = 0;
end

% Loesung des linearen GLS
T = A\F;

% Plot der Loesung
Z = zeros(Nx,Ny);
velocity = zeros(Nx,Ny);
for i = 1:Nx
    for j = 1:Ny
        Z(j,i) = T((j-1)*Nx+i); % quando quiser mudar o x,y troque i e j no Z(i,j). No momento está trocado
        velocity(j,i) = vx((j-1)*Nx+i); % aqui também está trocado
    end
end  

figure
contourf(Z',10)
title({'2-D Laplace''s Gleichung'})
xlabel('x-Achse')
ylabel('y-Achse')

figure
contourf(velocity',10)
title({'Updraft velocity'})
xlabel('x-Achse')
ylabel('y-Achse')

figure
surf(nodes_y,nodes_x,Z','EdgeColor','none');       
shading interp
title({'Moisture'})
xlabel('x-Achse (position,m)')
ylabel('y-Achse (height,m)')
zlabel('concentration [mol/m^3]')

