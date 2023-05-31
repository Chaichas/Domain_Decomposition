function [Ub, k, preSP_SP] = preCG_BDD(G, SP, bp, blockSD, blockAP, Ns, m)
    
    % preCG_BDD     - Preconditioned Conjugate gradient (Primal)
    %                       
    % Inputs:
    %   G           - Assembled rigid body modes
    %   SP          - Assembled primal Schur operator
    %   bp          - Assembled right-hand side of the primal interface problem
    %   blockSD     - Condensed dual Schur operators
    %   blockAP     - Condensed primal assembly operators
    %   Ns          - Number of subdomains
    %   m           - Maximum number of iterations
    %
    %
    % Outputs:
    %   Ub          - Computed displacement
    %   k           - Number of iterations to reach convergence
    %   preSP_SP    - multiplication of the assembled primal pre-conditioner with the assembled primal Schur

    %% set threshold for convergence criteria
    threshold = 1e-6;
    %% Compute the Neuwmann preconditioner
    Atild = inv(blockAP*blockAP')*blockAP; % it was assumed that blockM = identity matrix
    preSP = Atild * blockSD * Atild';

    preSP_SP = preSP*SP;

    %% Initialization
    G_invGtSPG_Gt = G*inv(G'*SP*G)*G';
    P = eye(Ns-1)-G_invGtSPG_Gt*SP;
    Ub = G_invGtSPG_Gt*bp;
    r = P'*bp;
    z = preSP*r;
    dtab = {};
    d = z;
    dtab{1} = d;
    norm_r0 = norm(r);
    ptab = {};
    for k = 0:m
        ptab{k+1} = P'*SP*d;
        %% Compute optimal step
        p = ptab{k+1};
        alpha = (r'*d)/(d'*p);
        %% Compute iterate
        Ub = Ub + alpha*d;
        %% Compute residual
        r = r - alpha*p;
        %% Check convergence criteria
        if norm(r)/norm_r0 < threshold
            break;
        end
        z = preSP*r;
        %% Reorthogonalisation
        beta = zeros(k+1);
        for j = 0:k
            beta(j+1) = -(z'*ptab{j+1})/(dtab{j+1}'*ptab{j+1});
        end
        %% Update search direction
        d = z;
        for j = 0:k
            d = d + beta(j+1)*dtab{k+1};
        end
        dtab{k+2} = d;
    end
end
