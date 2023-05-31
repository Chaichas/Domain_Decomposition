function [Ub, k] = CG_primal_interface(Ns,Ksub,Fsub,Ub,blockAP,m)

    % CG_primal_interface - Distributed Conjugate gradient method
    %                       
    % Inputs:
    %   Ns          - number of subdomains
    %   Ksub        - stiffness matrix on each subdomain
    %   Fsub        - load vector on each subdomain
    %   Ub          - initial displacement = 0
    %   blockAP     - concatenated primal Schur assembly operator
    %   m           - maximum number of iterations
    %
    % Outputs:
    %   Ub          - computed assembled displacment
    %   k           - number of iterations

    epsilon = 1e-6; % Threshold for convergence

    %% Compute local (concatenated) vector Ub(0)
    blockUb = blockAP'*Ub;
    last = 0;
    blockrb = zeros(size(blockUb));
    for s = 1:Ns      
        if (s == 1) || (s == Ns)
            Nb = 1;
        else
            Nb = 2;
        end
        % Solve Dirichlet problem
        Ui = inv(Ksub{s}(1:end-Nb,1:end-Nb))...
            *(Fsub{s}(1:end-Nb,1)-...
              Ksub{s}(1:end-Nb,end-Nb+1:end)*blockUb(last+1:last+Nb,1));
        % Compute local residual
        blockrb(last+1:last+Nb,1) = -Ksub{s}(end-Nb+1:end,1:end-Nb)*Ui-...
                                  Ksub{s}(end-Nb+1:end,end-Nb+1:end)*blockUb(last+1:last+Nb,1)+...
                                  Fsub{s}(end-Nb+1:end,1);
        last = last + Nb;
    end
    %% Compute global residual
    rb = blockAP*blockrb;
    db = rb;
    dbtab = {db};
    Spdbtab = {};
    norm_rb0 = norm(rb);

    for k = 0:m
        %% Compute local (concatenated) direction vector
        blockdb = blockAP'*db;
        last = 0;
        blockSpdb = zeros(size(blockdb));
        for s =1:Ns
            if (s == 1) || (s == Ns)
                Nb = 1;
            else
                Nb = 2;
            end
            % Solve local Dirichlet problem
            di = -inv(Ksub{s}(1:end-Nb,1:end-Nb))*...
                 Ksub{s}(1:end-Nb,end-Nb+1:end)*blockdb(last+1:last+Nb,1);
            %% Compute local matrix-vector product
            blockSpdb(last+1:last+Nb,1) = Ksub{s}(end-Nb+1:end,1:end-Nb)*di+...
                                        Ksub{s}(end-Nb+1:end,end-Nb+1:end)*blockdb(last+1:last+Nb,1);
            last = last + Nb;
        end
        %% Compute global matrix-vector product
        Spdb = blockAP*blockSpdb;
        Spdbtab{end+1} = Spdb;
        %% Compute optimal step
        alpha = (rb'*db)/(db'*Spdb);
        %% Compute iterate
        Ub = Ub + alpha*db;
        %% Compute residual
        rb = rb - alpha*Spdb;
        %% Reorthogonalization (from algorithm 3)
        beta = zeros(k+1);
        for j = 0:k
            beta(j+1) = -(rb'*Spdbtab{j+1})/(dbtab{j+1}'*Spdbtab{j+1});
        end
        %% Update search direction (from algorithm 3)
        db = rb;
        for j = 0:k
            db = db+beta(j+1)*dbtab{j+1};
        end
        dbtab{end+1} = db;
        %% Convergence criteria
        if norm(rb)/norm_rb0 < epsilon
            break
        end
    end
end