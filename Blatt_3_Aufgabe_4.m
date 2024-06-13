clear variables;
close all;

% Intervalle
Intx = [0 2];
Inty = [0 1];

% Dirichlet und Neumann Funktionen
gD = @(x,y) y*(1-y);
gN = @(x,y) -1;

% Quellterm
f = @(x,y) 0; 

n = 2;
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

%%% 5-Punkte Stern %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Aufstellen der Steifigkeitsmatrix
Np = Nx*Ny;
e = ones(Np,1);
A = -spdiags([hx*hx*e hy*hy*e -(2*hy*hy*e+2*hx*hx*e) hy*hy*e hx*hx*e],...
            [-Nx,-1,0,1,Nx],Np,Np);

% debugging
C = full(A);
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
for j = 2:Ny-1 % Linke Kante
    A((j-1)*Nx+1,:) = 0;
    A((j-1)*Nx+1,(j-1)*Nx+1) = 1;
    F((j-1)*Nx+1) = gD(0,(j-1)*hy);
end

% Einsetzen der Neumann-Randbedingungen
for j=2:Ny-1
    C = full(A);
    A((j-1)*Nx+Nx,:) = 0;
    C = full(A);
    A((j-1)*Nx+Nx,(j-1)*Nx+Nx-1) = -1;
    C = full(A);
    A((j-1)*Nx+Nx,(j-1)*Nx+Nx) = 1;
    F((j-1)*Nx+Nx) = hx*gN(2,(j-1)*hy);
end

% Loesung des linearen GLS
T = A\F;

%Plot der Loesung
Z = zeros(Nx,Ny);
for i = 1:Nx
    for j = 1:Ny
        Z(i,j) = T((j-1)*Nx+i); % quando quiser mudar o x,y troque i e j no Z(i,j). No momento está como no original
    end
end  

figure
contourf(Z',10)
title({'2-D Laplace''s Gleichung'})
xlabel('x-Achse')
ylabel('y-Achse')

figure
surf(nodes_x,nodes_y,Z','EdgeColor','none');       
shading interp
title({'2-D Laplace''s Gleichung'})
xlabel('x-Achse')
ylabel('y-Achse')
zlabel('Temperatur T')