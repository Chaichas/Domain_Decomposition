function [K,U] = unsubstructured(E,S,L,Fd,N)

    % unsubstructured  - Finite Element resolution of the 1D bar
    %                       
    % Inputs:
    %   E           - Young 's modulus
    %   S           - Section of the bar
    %   L           - Length of the bar
    %   Fd          - Tensile force
    %   N           - number of elements within each subdomain
    %
    % Outputs:
    %   K           - Stiffness matrix
    %   U           - Displacement

    
    dof = N+1; % Number of nodes
    h = L/N; % Size of each element

    K = zeros(dof,dof); % Initialization of the global stiffness matrix
    U = zeros(dof,1); % Initialization of the global displacement vector
    F = zeros(dof,1); % Initialization of the global force vector

    % Global stiffness Matrix
    ke = (S*E/h)*[1 -1; -1 1]; % Element stiffness matrix (uniform mesh)
    
    for i = 1:N
        var = [i i+1];
        K(var,var) = K(var,var) + ke; % Assembled global stiffness matrix
    end

    % Global force matrix
    F(dof,1) = Fd;

    % Elemination approach (using prescribed displacement boundary condition u(x=0)=0)

    % Eliminate the first row/ column of the global stiffness matrix K
    K(1,:) = []; % Eliminate the first row of the matrix K
    K(:,1) = []; % Eliminate the first column of the matrix K
    
    % Eliminate the first row of the global displacement vector D
    U(1,:) = []; % Eliminate the first row of D
    
    % Eliminate the first row of the global force vector F
    F(1,:) = []; % Eliminate the first row of F

    % System resolution
    U = K\F; % Solving of the system: [K][D]=[F]
end