%% Initialization (Question 1.2)
close all;
clear all;

% Input parameters
S = 600e-6; % Surface of the bar (in m^2)
E = 210e+9; % Young s modulus (in Pa)
L = 1; % Length of the bar (in m)
Fd = 1000; % Tensile force (in N)

N = 3; % Number of elements
dof = N+1; % Number of nodes

h = L/N; % Size of each element

K = zeros(dof,dof); % Initialization of the global stiffness matrix
D = zeros(dof,1); % Initialization of the global displacement vector
F = zeros(dof,1); % Initialization of the global force vector

%% Global stiffness Matrix

ke = (S*E/h)*[1 -1; -1 1]; % Element stiffness matrix (uniform mesh)

for i = 1:N
    var = [i i+1];
    K(var,var) = K(var,var) + (S*E/h)*[1 -1; -1 1]; % Assembled global stiffness matrix
end

%% Global force matrix

F(dof,1) = Fd;

%% Elemination approach (using boundary condition u(x=0)=0

% Eliminate the first row/ column of the global stiffness matrix K
K1 = K; % Create a copy of matrix K
K1(1,:) = []; % Eliminate the first row of the matrix K
K1(:,1) = []; % Eliminate the first column of the matrix K

% Eliminate the first row of the global displacement vector D
D1 = D; % create a copy of vector D
D1(1,:) = []; % Eliminate the first row of D

% Eliminate the first row of the global force vector F
F1 = F; % Create a copy of vector F
F1(1,:) = []; % Eliminate the first row of F

%% System resolution


D1 = K1\F1; % Solving of the system: [K][D]=[F]
D(2:dof) = D1; % Global displacement vector D

%% Display
% Global stiffness matrix
fprintf("******************************************* \n");
fprintf("The global stiffness matrix [K]: \n");
fprintf("******************************************* \n");
disp(K);
% Global force vector
fprintf("******************************************* \n");
fprintf("The global force vector [F]: \n");
fprintf("******************************************* \n");
disp(F);
% Global displacement vector
fprintf("******************************************* \n");
fprintf("The global displacement vector [D]: \n");
fprintf("******************************************* \n");
disp(D);

%% Plotting

% Coordinates along the x-axis
x = zeros(dof,1);
for i = 2:dof
    x(i) = x(i-1)+h;
end

% Displacement
figure
plot(x,D,'black--o',LineWidth = 3,MarkerSize= 12);
title('Displacement in function of x')
xlabel('Coordinate x [m]')
ylabel('Displacement [m]')
%grid on;

%% Strain and stress calculation

% Element strain vector
strain = zeros(N,1);
for i=1:N
    strain(i,1) = (1/h)*[-1 1]*[D(i,1);D(i+1,1)];
end

% Element stress vector
stress = E*strain;
