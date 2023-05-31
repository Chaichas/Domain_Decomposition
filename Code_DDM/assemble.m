function [SP, G, BP, blockAP, blockSP, blockSD, blockBP, blocke, blockBd, blockAD, Sd, Bd] = assemble(SPsub, SDsub, BPsub, Rbsub, APsub, ADsub, Ns)

    % assemble      - Assemble for the primal and dual approaches
    %                       
    % Inputs:
    %   SPsub       - Primal Shcur operator for each subdomain
    %   SDsub       - Dual Schur operatir for each subdomain
    %   BPsub       - Right-hand side of the primal interface problem for each subdomain
    %   Rbsub       - Rigid body modes for each subdomain
    %   APsub       - Primal assembly operator for each subdomain
    %   ADsub       - Dual assembly operator for each subdomain
    %   Ns          - Number of subdomains
    %
    %
    % Outputs:
    %   SP          - Assembled primal Schur operator
    %   G           - Assembled rigid body modes
    %   BP          - Assembled right-hand side of the primal interface problem
    %   blockAP     - Condensed primal assembly operators
    %   blockSP     - Condensed primal Schur operators
    %   blockSD     - Condensed dual Schur operators
    %   blockBP     - Condensed right-hand side of the primal interface problem
    %   blocke      - Condensed right-hand side of the dual interface problem
    %   blockBd     - Condensed right-hand side of the dual interface problem
    %   blockAD     - Condensed dual assembly operators
    %   Sd          - Assembled dual Schur operator
    %   Bd          - Assembled right-hand side of the dual interface problem



    %%
    % Dimension of concatenated matrix
    Nblock = 2*(Ns-1);

    % Initialization
    blockAP = zeros(Ns-1,Nblock); % Dimension of concatenated A (primal)
    blockAD = zeros(Ns-1,Nblock); % Dimension of concatenated Abar (dual)
    blockBP = zeros(Nblock,1); % Dimension of concatenated force
    blockSP = zeros(Nblock,Nblock); % Dimension of concatenated Sp
    blockSD = zeros(Nblock,Nblock); % Dimension of concatenated Sd
    blockRb = zeros(Nblock,Ns-1); % Dimension of concatenated rigid body mode matrix
    blockSD = zeros(Nblock,Nblock); % Dimension of concatenated Sd
    blockAD = zeros(Ns-1,Nblock); % Dimension of concatenated Abar (dual)
    last = 0;
    for s = 1:Ns
        if (s == 1) || (s == Ns)
            Nb = 1;
        else
            Nb = 2;
        end
        blockSP(last+1:last+Nb,last+1:last+Nb) = SPsub{s};
        blockSD(last+1:last+Nb,last+1:last+Nb) = SDsub{s};
        blockAP(:,last+1:last+Nb) = APsub{s};
        blockAD(:,last+1:last+Nb) = ADsub{s};
        blockBP(last+1:last+Nb,1) = BPsub{s};
        if (s == 1)
            blockRb(1,:) = 0; % 1st row of blockRb
        else
            blockRb(last+1:last+Nb, s-1) = Rbsub{s};
        end
        blockSD(last+1:last+Nb,last+1:last+Nb) = SDsub{s};
        blockAD(:,last+1:last+Nb) = ADsub{s};
        last = last + Nb;
    end
    blocke = blockRb'*blockBP;
    SP = blockAP*blockSP*blockAP'; % Assembled primal Schur complement
    BP = blockAP*blockBP; % Assembled right-hand side
    G = blockAD*blockRb;
    blockBd = blockSD*blockBP;
    Sd = blockAD*blockSD*blockAD.';
    Bd = blockAD*blockBd;

