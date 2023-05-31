function plot_U(Uisub, Ub, mapGlobalNodes, interfaceNode, L, Ns, N)

    % plot_U           - plot the computed displacement
    % approach
    %                       
    % Inputs:
    %   Uisub          - Displacement at internal nodes
    %   Ub             - Computed displacement
    %   mapGlobalNodes - Gloab nodes
    %   interfaceNode  - Interface nodes
    %   L              - Length of the bar
    %   Ns             - Number of subdomains
    %   N              - Number of elements

    %%    
    U = zeros(Ns*N,1);
    % Filling interface nodes
    for b = 1:Ns-1
        U(interfaceNode(b),1) = Ub(b);
    end
    % Filling internal nodes
    for s = 1:Ns
        if (s == 1) || (s == Ns)
            Nb = 1;
        else
            Nb = 2;
        end
        U(mapGlobalNodes{s}(1:end-Nb,1)',1) = Uisub{s};
    end

    % Computing the position of each node
    x = zeros(Ns*N+1,1);
    dx = L/(Ns*N);
    for i = 1:Ns*N+1
        x(i,1) = (i-1)*dx;
    end
    U;
    % plot
    plot(x, U)
end