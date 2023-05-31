function [Ub] = primalSchur(SP, BP)
   
    % decompose_primal_dual - Decomposition into subdomains 
    %                       
    % Inputs:
    %   SP          - Assembled primal Schur complement
    %   BP          - Assembled Right-hand side

    %
    % Outputs:
    %   Ub          - computed assembled displacment

    %% Direct resolution for the primal problem

    Ub = SP\BP; % Assembled displacement
end