function [Rbsub, blockRb] = rigid_body_mode(SPsub, Ns)

    % rigid_body_mode - Rigid body modes for the primal problem
    %                       
    % Inputs:
    %   SPsub       - Primal Schur complement for each subdomain
    %   Ns          - number of subdomains

    %
    % Outputs:
    %   Rbsub       - Rigid body mode for each subdomain
    %   blockRb     - Concatenated rigid body modes

    % Dimension of concatenated matrix
    Nblock = 2*(Ns-1);
    blockRb = zeros(Nblock,Ns-1); % Dimension of concatenated rigid body mode matrix
    
    Rbsub ={};
    last = 0;
    for s = 1:Ns
        if (s == 1) || (s == Ns)
            Nb = 1;
        else
            Nb = 2;
        end

        if (s == 1)
           Rbsub{s} = []; % last rigid body mode
           blockRb(1,:) = 0; % 1st row of blockRb
        elseif (s == Ns)
            Rbsub{s} = 1; % last rigid body mode
            blockRb(last+1:last+Nb, s-1) = Rbsub{s};
        else
           Rbsub{s} = [1;1];
           blockRb(last+1:last+Nb, s-1) = Rbsub{s};
        end
        last = last + Nb;
    end

% Rigid body modes using ker of Sp
%last = 0;
%     for s = 1:Ns
%         if (s == 1) || (s == Ns)
%             Nb = 1;
%         else
%             Nb = 2;
%         end
%         Rbsub{s} = null(SPsub{s});
%         if (s == 1)
%            blockRb(1,:) = 0; % 1st row of blockRb
%         else
%            blockRb(last+1:last+Nb, s-1) = Rbsub{s};
%         end
%         last = last + Nb;
%     end


end
