function [mapGlobalNodes, interfaceNode, APsub, ADsub, Ksub, SPsub, BPsub, Fsub, SDsub] = decompose_primal_dual(E,S,L,Fd,Ns,N)

    % decompose_primal_dual - Decomposition into subdomains (for primal/
    % dual)
    %                       
    % Inputs:
    %   E           - Young 's modulus
    %   S           - Section of the bar
    %   L           - Length of the bar
    %   Fd          - Tensile force
    %   Ns          - number of subdomains
    %   N           - number of elements within each subdomain

    %
    % Outputs:
    %   APsub       - Primal concatenated dual assembly operator
    %   ADsub       - Dual concatenated dual assembly operator
    %   Ksub        - stiffness matrix on each subdomain
    %   SPsub       - Primal Schur complement for each subdomain
    %   Fsub        - load vector on each subdomain
    %   BPsub       - Right-hand side vector
    %   SDsub       - Dual Schur complement for each subdomain

    H = L/Ns; % Length of the substruucture
    DOFsub = N+1; % Number of nodes in each substructure
    h = H/N; % Size of each element

    mapGlobalNodes = {}; % Converts local to global numbering
    globalNode = 1; % global node index
    interfaceNode = []; % Interface node
    % Definition of local and global numbering and primal assembly matrices
    for s = 1:Ns
        mapGlobalNodes{s} = zeros(DOFsub,1);
        for localNode = 1:DOFsub
            mapGlobalNodes{s}(localNode,1) = globalNode;
            globalNode = globalNode+1;
        end
        globalNode = globalNode-1;
        interfaceNode(s) = globalNode;
    end
    % Definition of the assembly operators for each substructure
    APsub = {}; % Primal assembly matrix for each substructure
    ADsub = {}; % DUal assembly matrix for each substructure
    APsub{1} = zeros(Ns-1,1);
    ADsub{1} = zeros(Ns-1,1);
    APsub{1}(1,1) = 1;
    ADsub{1}(1,1) = 1;
    for s = 2:Ns-1
        APsub{s} = zeros(Ns-1,2);
        APsub{s}(s-1,1) = 1;
        APsub{s}(s,2) = 1;
        ADsub{s} = zeros(Ns-1,2);
        ADsub{s}(s-1,1) = -1;
        ADsub{s}(s,2) = 1;
    end
    APsub{Ns} = zeros(Ns-1,1);
    APsub{Ns}(end,1) = 1;
    ADsub{Ns} = zeros(Ns-1,1);
    ADsub{Ns}(end,1) = -1;
    
    Fsub = {}; % Force vector for each substructure
    Usub = {}; % Displacement vector for each substructure
    Ksub = {}; % Stifness matrix for each substructure

    ke = (S*E/h)*[1 -1; -1 1]; % Element stiffness matrix (uniform mesh)

    % Definition of stiffness matrix and forces for each substructure
    for s = 1:Ns
        Ksub{s} = zeros(DOFsub,DOFsub);
        Fsub{s} = zeros(DOFsub,1);
        Usub{s} = zeros(DOFsub,1);
        for localElm = 1:N
            var = [localElm localElm+1];
            Ksub{s}(var,var) = Ksub{s}(var,var) + ke;
        end
    end
    Fsub{Ns}(DOFsub,1) = Fd;

    %% Substructure 1 (Prescribed Displacement Boundary condition)
    % Elimination approach (using prescribed displacement boundary condition u(x=0)=0)

    % Eliminate the first row/ column of the global stiffness matrix K
    Ksub{1}(:,1) = []; % Eliminate the first row of the matrix K
    Ksub{1}(1,:) = []; % Eliminate the first column of the matrix K
    
    % Eliminate the first row of the global displacement vector U
    Usub{1}(1,:) = []; % Eliminate the first row of D
    
    % Eliminate the first row of the global force vector F
    Fsub{1}(1,:) = []; % Eliminate the first row of F
    
    mapGlobalNodes{1}(1,:) = []; % Eliminate 1st node from local and global nodes

    %% Reorganise the rows of the previous matrix and vectors in internal and interface nodes

    % The first subdomain is already organized
    % For the last subdomain, we should move the first node to the last position
    % For the other subdomains, we need to move the first node to the N-1 position
    
    for s = 2:Ns
        if s == Ns
            last = DOFsub;
        else
            last = DOFsub-1;
        end
        % Swap rows of the stiffnes matrix
        tmp = Ksub{s}(:,1);
        for col = 2:last
            Ksub{s}(:,col-1) = Ksub{s}(:,col);
        end
        Ksub{s}(:,last) = tmp;
        % Swap columns of the stiffnes matrix
        tmp = Ksub{s}(1,:);
        for row = 2:last
            Ksub{s}(row-1,:) = Ksub{s}(row,:);
        end
        Ksub{s}(last,:) = tmp;
        % Swap rows of the displacement vector
        tmp = Usub{s}(1,1);
        for row = 2:last
            Usub{s}(row-1,1) = Usub{s}(row,1);
        end
        Usub{s}(last,1) = tmp;
        % Swap rows of the force vector
        tmp = Fsub{s}(1,1);
        for row = 2:last
            Fsub{s}(row-1,1) = Fsub{s}(row,1);
        end
        Fsub{s}(last,1) = tmp;
        % Also swap row of mapGlobalNode
        tmp = mapGlobalNodes{s}(1,1);
        for row = 2:last
            mapGlobalNodes{s}(row-1,1) = mapGlobalNodes{s}(row,1);
        end
        mapGlobalNodes{s}(last,1) = tmp;
    end

    %% Compute the Primal Schur complement, rigid-body modes and primal right-hand side of each subdomain
    SPsub = {};
    Rbsub = {};
    for s = 1:Ns
        if (s == 1) || (s == Ns)
            Nb = 1;
        else
            Nb = 2;
        end
        KbiinvKii = Ksub{s}(end-Nb+1:end,1:end-Nb)*inv(Ksub{s}(1:end-Nb,1:end-Nb));
        SPsub{s} = Ksub{s}(end-Nb+1:end,end-Nb+1:end) ...
                      - KbiinvKii*Ksub{s}(1:end-Nb,end-Nb+1:end);
        %Rbsub{s} = null(SPsub{s});
        BPsub{s} = Fsub{s}(end-Nb+1:end,1) ...
                      - KbiinvKii*Fsub{s}(1:end-Nb,1);
    end

    %% Compute the Dual Schur complement
    SDsub = {};
    for s = 1:Ns
        SDsub{s} = pinv(SPsub{s}); % Pseudo-inverse
    end

end