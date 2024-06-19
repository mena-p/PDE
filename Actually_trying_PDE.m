clear variables;
close all;
%% README 

% Hi, Nils. This is the code used to generate the results in the write-up.
% You can alter the simulation and physical parameters in their respective
% sections. You don't need to change anything else in the code.

% The set-up section creates the grid and the vector with velocities in
% the lexicographic order, which we need to create the matrix implementing
% the five difference star in the section after that. The matrix is altered
% in the section "boundary conditions"  to account for... boundary 
% conditions. The last section solves the equation ad plots everything.

%% Pre-defined cases

% I included some presets to make things easier. To reproduce the results
% in the write-up, choose a preset from 1 to 9. Presets 10 and 11 are for
% validating the program by solving the advection difusion equation on a 1
% by 1 grid.
% 
% Choose 0 if you wish to play around with the parameters yourself.

PRESET = 3;

%% Simulation parameters

% Domain size (meters)
intx = [20 300]; % this is the height
inty = [-50 50]; % this is the horizontal position

% Set the amount of grid points (change only the n)
n = 8;
nx = 2^n;
ny = 2^n;
%% Physical parameters

% Source therm
f = @(x,y) 0; % zero since there is no chemical reaction

% Diffusivity of water vapor in air
D = 0.000025; % m^2/s

% Boundary conditions
% note: any function of x, y that returns a scalar can be used as a 
% neumann condition. Other possible choices are BC_normal_FD(x,y) or 
% BC_discontinous_FD(x,y), but they yield crappy results.

dirichlet = @(x,y) 0.5; % Humidity concentration on the left and right 
                        % edges of the domain, mol/m^3

% a value 100000x smaller than the one I derived in the write-up yields
% physical results :/

neumann_top = @(x,y) -0.00413; % Humidity concentration gradient (flux) on 
                               % top of domain, (mol/m^3)/m
neumann_sur = @(x,y) -0.00413; % Humidity concentration gradient (flux) at 
                               % the surface, (mol/m^3)/m

% Updraft velocity function
% note: you can choose any function of x, y that returns a scalar, but 
% keep in mind that y=location at surface and x=height in this 
% implementation.

velocity = @(height,position) allen(position,height,4,1000,250000,0);

%% Parameters for the pre-defined cases

% Domain size (meters)
intx = [20 500]; % this is the height
inty = [-50 50]; % this is the horizontal position

% Set the amount of grid points
n = 8;
nx = 2^n;
ny = 2^n;

% Source therm
f = @(x,y) 0;

% Diffusivity of water vapor in air
D = 0.000025; % m^2/s

if PRESET == 1

    % Boundary conditions
    dirichlet = @(x,y) 0.5;
    neumann_top = @(x,y) 0;
    neumann_sur = @(x,y) 0;
   
    velocity = @(height,position) allen(position,height,4,1000,250000,0);

elseif PRESET == 2

    % Boundary conditions
    dirichlet = @(x,y) 0.5;
    neumann_top = @(x,y) -413;
    neumann_sur = @(x,y) -413;
    
    velocity = @(height,position) allen(position,height,4,1000,250000,0);

elseif PRESET == 3
    
    % Boundary conditions
    dirichlet = @(x,y) 0.5;
    neumann_top = @(x,y) BC_normal_FD(x,y);
    neumann_sur = @(x,y) BC_normal_FD(x,y);
    
    velocity = @(height,position) allen(position,height,4,1000,250000,0);

elseif PRESET == 4
    
    % Boundary conditions
    dirichlet = @(x,y) 0.5;
    neumann_top = @(x,y) -413;
    neumann_sur = @(x,y) 0;
    
    velocity = @(height,position) allen(position,height,4,1000,250000,0);

elseif PRESET == 5
    
    % Boundary conditions
    dirichlet = @(x,y) 0.5;
    neumann_top = @(x,y) BC_normal_FD(x,y);
    neumann_sur = @(x,y) 0;
    
    velocity = @(height,position) allen(position,height,4,1000,250000,0);

elseif PRESET == 6
    
    % Boundary conditions
    dirichlet = @(x,y) 0.5;
    neumann_top = @(x,y) 0;
    neumann_sur = @(x,y) -413;
    
    velocity = @(height,position) allen(position,height,4,1000,250000,0);

elseif PRESET == 7
    
    % Boundary conditions
    dirichlet = @(x,y) 0.5;
    neumann_top = @(x,y) BC_normal_FD(x,y);
    neumann_sur = @(x,y) -413;
    
    velocity = @(height,position) allen(position,height,4,1000,250000,0);

elseif PRESET == 8
    
    % Boundary conditions
    dirichlet = @(x,y) 0.5;
    neumann_top = @(x,y) 0;
    neumann_sur = @(x,y) BC_normal_FD(x,y);
    
    velocity = @(height,position) allen(position,height,4,1000,250000,0);

elseif PRESET == 9

    % Boundary conditions
    dirichlet = @(x,y) 0.5;
    neumann_top = @(x,y) -413;
    neumann_sur = @(x,y) BC_normal_FD(x,y);
    
    velocity = @(height,position) allen(position,height,4,1000,250000,0);

elseif PRESET == 10

    % Domain size (meters)
    intx = [0 1];
    inty = [-0.5 0.5];
    
    % Diffusivity
    D = 0.25; % m^2/s
    
    % Boundary conditions
    dirichlet = @(x,y) 0.5;
    neumann_top = @(x,y) -0.5;
    neumann_sur = @(x,y) -0.5;
    
    velocity = @(height,position) 0;

elseif PRESET == 11

    % Domain size (meters)
    intx = [0 1];
    inty = [-0.5 0.5];

    % Diffusivity
    D = 0.25; % m^2/s
    
    % Boundary conditions
    dirichlet = @(x,y) 0.5;
    neumann_top = @(x,y) -0.5;
    neumann_sur = @(x,y) -0.5;
    
    velocity = @(height,position) 1;
end
%% Set up

% Grid spacing
hx = (intx(2)-intx(1))/nx;
hy = (inty(2)-inty(1))/ny;

% Creating the vector of grid nodes and including points in the edge
Nx = nx+1;
Ny = ny+1;
nodes_x = linspace(intx(1),intx(2),Nx);
nodes_y = linspace(inty(1),inty(2),Ny);


% Create vector with updraft velocities (lexicographic order)
vx = zeros(Nx*Ny,1);
for i = 1:Nx
    for j = 1:Ny
        % coordinate where velocity is calculated
        x = nodes_x(i);
        y = nodes_y(j);

        % lexicographic index
        l = (j-1)*Nx+i;
        vx(l) = 1*velocity(x,y); % y is location, x is height in our PDE!
    end
end

% correct lexicographic order so that it fits the -1 diagonal of the
% matrix A (by removing the first element of the array and shifting the
% values one index right)
vx(1:Nx*Ny-1) = vx(2:Nx*Ny);

% remove NaN (some velocity functions may return NaN)
vx(isnan(vx)) = 0;

%% 5-Point difference star

% Set up the finite difference matrix
Np = Nx*Ny; % number of grid points
e = ones(Np,1); % unity vector

A =  spdiags([D*hx*hx*e D*hy*hy*e+hy*hy*hx.*vx ...
    -2*D*(hy*hy+hx*hx)*e-hx*hy*hy*vx D*hy*hy*e D*hx*hx*e],...
            [-Nx,-1,0,1,Nx],Np,Np);

% Set up of the RHS of the equation (just zero in our case, changing the 
% source function changes this)
F = zeros(Np,1);
for j = 1:Ny
    for i = 1:Nx
        F((j-1)*Nx+i) = hx*hx*hy*hy*f((i-1)*hx,(j-1)*hy); 
    end
end

%% Boundary conditions

% Imposing Neumann boundary conditions
start_x = nodes_x(1,1);
start_y = nodes_y(1,1);
for j=2:Ny-1 % Right edge of domain (top of thermal)
    A((j-1)*Nx+Nx,:) = 0;
    % Backwards difference to approximate derivative
    A((j-1)*Nx+Nx,(j-1)*Nx+Nx-1) = -1;
    A((j-1)*Nx+Nx,(j-1)*Nx+Nx) = 1;
    F((j-1)*Nx+Nx) = hx*neumann_top(start_x+(Nx-1)*hx, start_y+(j-1)*hy);
end
for j=2:Ny-1 % Left edge of domain (surface)
    A((j-1)*Nx+1,:) = 0;
    % Forwards difference to approximate derivative
    A((j-1)*Nx+1,(j-1)*Nx+1) = -1;
    A((j-1)*Nx+1,(j-1)*Nx+1+1) = 1;
    F((j-1)*Nx+1) = hx*neumann_sur(start_x+(1-1)*hx, start_y+(j-1)*hy);
end

% Imposing Dirichlet boundary conditions
for i = 1:Nx % Lower edge (right side of thermal)
    A(i,:) = 0;
    A(i,i) = 1;
    F(i) = 0;
end
for i = 1:Nx % Upper edge (left side of thermal)
    A((Ny-1)*Nx+i,:) = 0;
    A((Ny-1)*Nx+i,(Ny-1)*Nx+i) = 1;
    F((Ny-1)*Nx+i) = 0;
end
%% Solution and plotting

% Solving the system of linear equations
T = A\F;

% Get result out of vector in lexicographic order and into the matrix C
C = zeros(Nx,Ny); % humidity concentration
vel = zeros(Nx,Ny); % updraft velocity
for i = 1:Nx
    for j = 1:Ny
        C(j,i) = T((j-1)*Nx+i); % j and i are flipped to undo the choice 
                                % of y = location and x = height
        vel(j,i) = vx((j-1)*Nx+i); % here as well
    end
end  

% Plot everything

flux = neumann_sur(0,nodes_y);
figure
plot(nodes_y,flux)
title('Boundary condition on surface')
xlabel('Position')
ylabel('dc/dx')

% Humidity countour plot
figure
contourf(nodes_y,nodes_x,C',10)
title('Water vapor')
xlabel('position, m')
ylabel('height, m')
colorbar

% Velocity countour plot
figure
contourf(nodes_y,nodes_x,vel',10)
axis equal
title('Updraft velocity')
xlabel('position, m')
ylabel('height, m')

% Humidity surface plot
figure
surf(nodes_y,nodes_x,C','EdgeColor','none');       
shading interp
title('Water vapor')
xlabel('position, m')
ylabel('height, m')
zlabel('concentration [mol/m^3]')
