function Ub = dualSchur(SDsub, BPsub, Rbsub, blocke, blockAD, Ns, Sd, Bd, G)

    % dualSchur     - Resolution of the Schur interface problem using a
    % direct method
    %                       
    % Inputs:
    %   SDsub       - Dual Schur operatir for each subdomain
    %   BPsub       - ight-hand side of the primal interface problem for each subdomain
    %   Rbsub       - Primal Schur complement for each subdomain
    %   blocke      - Condensed right-hand side of the dual interface problem
    %   blockAD     - Condensed dual assembly operators
    %   Ns          - Number of subdomains
    %   Sd          - Assembled dual Schur operator
    %   Bd          - Assembled right-hand side of the dual interface problem
    %   G           - Assembled rigid body modes

    %
    % Outputs:
    %   Ub          - Computed displacement

    %% Assemble the global primal Schpur complement

    % Dimension of concatenated matrix
    Nblock = 2*(Ns-1);

    L = zeros(Nblock,Nblock);
    L(1:Ns-1,1:Ns-1) = Sd;
    L(1:Ns-1,Ns:end) = G;
    L(Ns:end,1:Ns-1) = G.';

    RHS = zeros(Nblock,1);
    RHS(1:Ns-1,1) = -Bd;
    RHS(Ns:end,1) = -blocke;

    Sol = L\RHS;
    stress = Sol(1:Ns-1,1);
    alpha = Sol(Ns:end,1);

    block_stress = blockAD'*stress;

    stress_sub = {};
    alpha_sub = {};
    Ubsub = {};
    last = 0;
    for s = 1:Ns
        if (s == 1) || (s == Ns)
            Nb = 1;
        else
            Nb = 2;
        end

        stress_sub{s} = block_stress(last+1:last+Nb,1);

        if s == 1
            alpha_sub{s} = 0;
            Ubsub{s} = SDsub{s}*(BPsub{s}+stress_sub{s});
        else
            alpha_sub{s} = alpha(s-1);
            Ubsub{s} = SDsub{s}*(BPsub{s}+stress_sub{s})+Rbsub{s}*alpha_sub{s};
        end

        last = last +1;
    end
    
    Ub = zeros(Ns-1,1);
    Ub(1,1) = Ubsub{1};
    for s = 2:Ns-1
        Ub(s,1) = Ubsub{s}(end,1);
    end
end
