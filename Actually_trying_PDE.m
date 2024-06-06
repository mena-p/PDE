clear variables;
close all;

%% Simulation parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Domain size (meters)
intx = [200 500];  % this is the height
inty = [-100 100]; % this is the horizontal position

% Set the amount of grid points (change only the n)
n = 9;
nx = 2^n;
ny = 2^n;
%% Physical parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Source therm
f = @(x,y) 0; % zero since there is no chemical reaction

% Diffusivity of water vapor in air
D = 0.000025; % m^2/s

% Boundary conditions
dirichlet = @(x,y) 0; % Humidity concentration on the left and right edges of the domain, mol/m^3
neumann = @(x,y) 400; % Humidity concentration gradient (flux) on top and bottom of domain, (mol/m^3)/m
Boundary_condition = @(x,y) 400;

% Choosen updraft velocity function
% note: you can choose any function of x,y that returns a scalar. However,
% keep in mind that y = location at surface and x = height in this implementation.

velocity = @(height,position) allen(position,height,4,1000,250000,0);

%% Set up %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%% 5-Point difference star %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up the finite difference matrix
Np = Nx*Ny; % number of grid points
e = ones(Np,1); % unity vector

A =  spdiags([D*hx*hx*e D*hy*hy*e+hy*hy*hx.*vx -2*D*(hy*hy+hx*hx)*e-hx*hy*hy*vx D*hy*hy*e D*hx*hx*e],...
            [-Nx,-1,0,1,Nx],Np,Np);

%%% Advective part (removed and included directly in the code above)
% a = 1;
% b = 12;
% B = spdiags([-0.5*b*hy*hx*hx*e -0.5*a*hx*hy*hy*e 0.5*a*hx*hy*hy*e 0.5*b*hy*hx*hx*e],...
%             [-Nx,-1,1,Nx],Np,Np);
% A = A + B;

% Set up of the RHS of the equation (just zero in our case, changing the 
% source function changes this)
F = zeros(Np,1);
for j = 1:Ny
    for i = 1:Nx
        F((j-1)*Nx+i) = hx*hx*hy*hy*f((i-1)*hx,(j-1)*hy); 
    end
end

%% Boundary conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Imposing Neumann boundary conditions
for j=2:Ny-1 % Right edge of domain (top of thermal)
    A((j-1)*Nx+Nx,:) = 0;
    % Backwards difference to approximate derivative
    A((j-1)*Nx+Nx,(j-1)*Nx+Nx-1) = -1;
    A((j-1)*Nx+Nx,(j-1)*Nx+Nx) = 1;
    F((j-1)*Nx+Nx) = -hx*neumann((Nx-1)*hx,(j-1)*hy);
end
for j=2:Ny-1 % Left edge of domain (surface)
    A((j-1)*Nx+1,:) = 0;
    % Forwards difference to approximate derivative
    A((j-1)*Nx+1,(j-1)*Nx+1) = -1;
    A((j-1)*Nx+1,(j-1)*Nx+1+1) = 1;
    F((j-1)*Nx+1) = -hx*neumann((1-1)*hx,(j-1)*hy);
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
        C(j,i) = T((j-1)*Nx+i); % j and i are flipped to undo the choice of y = location and x = height
        vel(j,i) = vx((j-1)*Nx+i); % here as well
    end
end  

% Plot everything

% Humidity countour plot
figure
contourf(nodes_y,nodes_x,C',10)
title({'2-D Laplace''s Gleichung'})
xlabel('x-Achse')
ylabel('y-Achse')

% Velocity countour plot
figure
contourf(nodes_y,nodes_x,vel',10)
title({'Updraft velocity'})
xlabel('x-Achse')
ylabel('y-Achse')

% Humidity surface plot
figure
surf(nodes_y,nodes_x,C','EdgeColor','none');       
shading interp
title({'Moisture'})
xlabel('x-Achse (position,m)')
ylabel('y-Achse (height,m)')
zlabel('concentration [mol/m^3]')

