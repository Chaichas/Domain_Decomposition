function [Ub, k] = FETI_PCPG(G, SDsub, SPsub, BPsub, RBsub, ADsub, blockBd, blocke, blockAD, blockSP, Ns, m)
    
    % FETI_PCPG     - Distributed FETI algorithm (Projected CG) for dual
    % approach
    %                       
    % Inputs:
    %   G           - Assembled rigid body modes
    %   SDsub       - Dual Schur operatir for each subdomain
    %   SPsub       - Primal Schur complement for each subdomain
    %   BPsub       - Right-hand side of the primal interface problem for each subdomain
    %   RBsub       - Rigid body modes for each subdomain
    %   ADsub       - Dual assembly operator for each subdomain
    %   blockBd     - Condensed right-hand side of the dual interface problem
    %   blocke      - Condensed right-hand side of the dual interface problem
    %   blockAD     - Condensed dual assembly operators
    %   blockSP     - Condensed primal Schur operators
    %   Ns          - Number of subdomains
    %   m           - Maximum number of iterations

    %
    % Outputs:
    %   Ub          - Computed displacement
    %   k           - Number of iterations to reach convergence

    %% Dimension of concatenated matrix
    Nblock = 2*(Ns-1);

    % set threshold for convergence criteria
    threshold = 1e-6;
    
    %% Compute the preconditioner for each subdomain
    preSD = zeros(Ns-1,Ns-1);
    %M = 0.5*diag(Ns-1); % It is assumed that the multiplicity is equal to 2
    %for s = 1:Ns
    %    preSD = preSD + M*ADsub{s}*SPsub{s}*ADsub{s}'*M;
    %end
    Atild = inv(blockAD*blockAD')*blockAD;
    preSD = Atild * blockSP * Atild';
    
    % Initialization

    %% Compute P and initial lambda
    invGtG = inv(G'*G);
    P = eye(Ns-1)-G*invGtG*G';
    lam = -G*invGtG*blocke; % Initial solution

    %% Compute initial residual (distributed)
    % Distribute the lam
    block_lam = blockAD'*lam;
    % Compute Sd*lam at the subdomain level
    block_bd_Sdlam = zeros(Nblock,1);
    last = 0;
    for s = 1:Ns
        if (s == 1) || (s == Ns)
            Nb = 1;
        else
            Nb = 2;
        end
        block_bd_Sdlam(last+1:last+Nb,1) = -blockBd(last+1:last+Nb,1) - SDsub{s}*block_lam(last+1:last+Nb,1);
        last = last + Nb;
    end
    % Assemble -bd-Sd*lam and multiply with P^T
    r = P'*blockAD*block_bd_Sdlam;
    norm_r0 = norm(r);
    k=0;
    if (norm_r0 > threshold)
        %% Compute z
        z = P*preSD*r;
        % Set the initial search direction
        d = z;
        %% start loop
        dtab = {};
        ptab = {};
        dtab{1} = d;
        for k = 0:m
            %% Set p (distributed)
            % Distribute d
            block_d = blockAD'*d;
            % Compute Sd*d for each subdomain
            block_Sdd = zeros(Nblock,1);
            last = 0;
            for s = 1:Ns
                if (s == 1) || (s == Ns)
                    Nb = 1;
                else
                    Nb = 2;
                end
                block_Sdd(last+1:last+Nb,1) = SDsub{s}*block_d(last+1:last+Nb,1);
                last = last + Nb;
            end
            % Assemble Sd*d and compute p
            p = P'*blockAD*block_Sdd;
            ptab{k+1} = p;
            %% Compute optimal step
            alpha = (r'*d)/(d'*p);
            %% Compute iterate
            lam = lam + alpha*d;
            %% compute residual
            r = r - alpha*p;
            %% Check convergence
            if norm(r)/norm_r0 < threshold
                break;
            end
            %% Set z
            z = P*preSD*r;
            %% Reorthogonalization
            beta = zeros(k+1);
            for j = 0:k
                beta(j+1) = -(z'*ptab{j+1})/(dtab{j+1}'*ptab{j+1});
            end
            %% Update search direction
            d = z;
            for j = 0:k
                d = d + beta(j+1)*dtab{j+1};
            end
            dtab{k+2} = d;
        end % end loop
    end
    %% Post-processing
    
    %% Compute block_alpha (distributed)
    % Distribute lambda
    block_lam = blockAD'*lam;
    % Compute Sd*lam for each subdomain
    block_bd_Sdlam = zeros(Nblock,1);
    last = 0;
    for s = 1:Ns
        if (s == 1) || (s == Ns)
            Nb = 1;
        else
            Nb = 2;
        end
        block_bd_Sdlam(last+1:last+Nb,1) = -blockBd(last+1:last+Nb,1) - SDsub{s}*block_lam(last+1:last+Nb,1);
        last = last + Nb;
    end
    block_alphab = invGtG*G'*blockAD*block_bd_Sdlam;

    %% Compute Ub
    Ubsub = {};
    block_lam = blockAD'*lam;
    last = 0;
    for s = 1:Ns
        if (s == 1) || (s == Ns)
            Nb = 1;
        else
            Nb = 2;
        end

        lam_s = block_lam(last+1:last+Nb,1);

        if s == 1
            Ubsub{s} = SDsub{s}*(BPsub{s}+lam_s);
        else
            Ubsub{s} = SDsub{s}*(BPsub{s}+lam_s)+RBsub{s}*block_alphab(s-1);
        end
        last = last +1;
    end

    Ub = zeros(Ns-1,1);
    Ub(1,1) = Ubsub{1};
    for s = 2:Ns-1
        Ub(s,1) = Ubsub{s}(end,1);
    end
end
