close all;
clear all;

%% Input parameters

S = 600e-6; % Surface of the bar (in m^2)
E = 210e+9; % Young s modulus (in Pa)
L = 1; % Length of the bar (in m)
Fd = 1000; % Tensile force (in N)

Ns = 3; % Number of substructures
N = 2; % Number of elements in each substructure

%% (Q1.3) : Primal and Dual subdomain decomposition
[mapGlobalNodes, interfaceNode, APsub, ADsub, Ksub, SPsub, BPsub, Fsub, SDsub] = decompose_primal_dual(E,S,L,Fd,Ns,N);
[Rbsub, blockRb] = rigid_body_mode(SPsub, Ns); % rigid body modes
[SP, G, BP, blockAP, blockSP, blockSD, blockBP, blocke, blockBd, blockAD, Sd, Bd] = assemble(SPsub, SDsub, BPsub, Rbsub, APsub, ADsub, Ns);
%% (Q2.1) : Direct resolution of the primal Schur interface problem
[Ub_primal_direct] = primalSchur(SP, BP);

%% (Q2.2) : Evaluation of the conditioning

% Conditioning of Sp with respect to h/H
% The sensitivity to h/H is equivalent to the sensitivity to 1/N
% We will fix Ns = 5 (H fixed) and vary N from 2 to 20 by step of 2

% Ns = 10;
% condSP = []; % Contditioning of assembled Schur matrix
% condK = []; % Conditioning of stiffness matrix (unsubstructured)
% Ntab = 2:2:50;
% rNtab = [];
% for N = Ntab
%     [mapGlobalNodes, interfaceNode, APsub, ADsub, Ksub, SPsub, BPsub, Fsub, SDsub] = decompose_primal_dual(E,S,L,Fd,Ns,N);
%     [Rbsub, blockRb] = rigid_body_mode(SPsub, Ns); % rigid body modes
%     [SP, G, BP, blockAP, blockSP, blockSD, blockBP, blocke, blockBd, blockAD, Sd, Bd] = assemble(SPsub, SDsub, BPsub, Rbsub, APsub, ADsub, Ns);
%     [Ub_primal_direct] = primalSchur(SP, BP);
%     [K,U] = unsubstructured(E,S,L,Fd,N*Ns);
%     lambdaSP = abs(eig(SP));
%     lambdaK = abs(eig(K));
%     condSP(end+1) = max(lambdaSP)/min(lambdaSP);
%     condK(end+1) = max(lambdaK)/min(lambdaK);
%     rNtab(end+1) = 1/N;
% end
% figure
% plot(rNtab,condSP, 'red-o',rNtab, condK, 'blue-o', 'LineWidth', 2,MarkerSize= 8);
% title('Evolution of the conditioning of assembled Sp and K in function of the parameter (h/H)')
% xlabel('Parameter (h/H = 1/N) with Ns fixed', 'FontWeight', 'bold')
% ylabel('Conditioning', 'FontWeight', 'bold')
% legend('condSP', 'condK', 'FontWeight', 'bold', 'FontSize', 14)
% 
%% (Q2.3) : Resolution with the distributed CG (iterative)
Ub_primal_iter = zeros(Ns-1,1);
m = 100;
% Start timer
tic;
[Ub_primal_iter,k_primalGC] = CG_primal_interface(Ns,Ksub,Fsub,Ub_primal_iter,blockAP,m);
% Stop timer and display elapsed time
elapsed_time_CG_primal = toc;
fprintf('Elapsed time of the distributed CG algorithm: %f seconds\n', elapsed_time_CG_primal);
%fprintf('Number of iterations: %d \n', k_primalGC);

%% Plotting of the displacement in the primal approach (possible for all approaches)
Uisub = compute_Ui(Ksub, Fsub, Ub_primal_iter, blockAP, Ns);
% plot_U(Uisub, Ub_primal_iter, mapGlobalNodes, interfaceNode, L, Ns, N);
% title('Displacement in function of the x-coordinate');
% xlabel('x-coordinate (m)');
% ylabel('Displacement')

%% (Q2.4) : Numerical Scalability of distributed CG

% m = 100;
% iter = []; % Number of iterations
% time = []; % Computational time
% N = 3; % Number of elements (fixed)
% Ntab1= 2:2:200; % Number of subdomains
% H = [];
% rNtab = [];
% for Ns = Ntab1
%     [mapGlobalNodes, interfaceNode, APsub, ADsub, Ksub, SPsub, BPsub, Fsub, SDsub] = decompose_primal_dual(E,S,L,Fd,Ns,N);
%     [Rbsub, blockRb] = rigid_body_mode(SPsub, Ns); % rigid body modes
%     [SP, G, BP, blockAP, blockSP, blockSD, blockBP, blocke, blockBd,blockAD, Sd, Bd] = assemble(SPsub, SDsub, BPsub, Rbsub, APsub, ADsub, Ns);
%     Ub_iter = zeros(Ns-1,1);
%     tic;
%     [Ub_iter,k_primalGC] = CG_primal_interface(Ns,Ksub,Fsub,Ub_iter,blockAP,m);
%     t = toc;
%     iter(end+1) = k_primalGC;
%     time(end+1) = t;
%     rNtab(end+1) = Ns;
%     H(end+1) = L / Ns; % length of subdomains
% end
% % %Create figure 
% %figure
% % 
% % % %Define subplot 1 
% % subplot(1,2,1);
% % plot(rNtab,iter, 'blue-o', 'LineWidth',2,MarkerSize= 8);
% % title('Evolution of the number of iterations in function of number of subdomains', 'FontWeight', 'bold') ;
% % xlabel('Number of subdomains', 'FontWeight', 'bold') ;
% % ylabel('Number of iterations', 'FontWeight', 'bold');
% % 
% % % Define subplot 2 
% % subplot(1,2,2);
% % plot(rNtab,time, 'blue-o', 'LineWidth',2,MarkerSize= 8);
% % title('Evolution of the Computational time in function of number of subdomains', 'FontWeight', 'bold');
% % xlabel('Number of subdomains', 'FontWeight', 'bold');
% % ylabel('Computational time (s)', 'FontWeight', 'bold');
% 
% %Create figure 2
% figure
% 
% %Define subplot 1
% plot(H,iter, 'red-o', 'LineWidth', 2,MarkerSize= 8)
% title('Evolution of the number of iterations in function of H', 'FontWeight', 'bold')
% xlabel('Length of subdomain H', 'FontWeight', 'bold')
% ylabel('Number of iterations', 'FontWeight', 'bold')

%% (Q2.5) : Resolution with preconditioned CG
m = 100; % Maximal number of iterations
tic;
[Ub_primal_iter_pre, k_preGC, preSP_SP] = preCG_BDD(G, SP, BP, blockSD, blockAP, Ns, m);
elapsed_time_pre = toc;
fprintf('Elapsed time of the primal preconditioner algorithm: %f seconds\n', elapsed_time_pre);


%% (Q2.5) : Numerical Scalability of preconditioned CG

% m = 100;
% iter = []; % Number of iterations
% time = []; % Computational time
% N = 3; % Number of elements (fixed)
% Ntab1= 2:2:200; % Number of subdomains
% H = [];
% rNtab = [];
% for Ns = Ntab1
%     [mapGlobalNodes, interfaceNode, APsub, ADsub, Ksub, SPsub, BPsub, Fsub, SDsub] = decompose_primal_dual(E,S,L,Fd,Ns,N);
%     [Rbsub, blockRb] = rigid_body_mode(SPsub, Ns); % rigid body modes
%     [SP, G, BP, blockAP, blockSP, blockSD, blockBP, blocke, blockBd,blockAD, Sd, Bd] = assemble(SPsub, SDsub, BPsub, Rbsub, APsub, ADsub, Ns);
%     Ub_iter = zeros(Ns-1,1);
%     tic;
%     [Ub, k_preGC, preSP_SP] = preCG_BDD(G, SP, BP, blockSD, blockAP, Ns, m);
%     t = toc;
%     iter(end+1) = k_preGC;
%     time(end+1) = t;
%     rNtab(end+1) = Ns;
%     H(end+1) = L / Ns; % length of subdomains
% end
% % Create figure 
% figure
% 
% % number of iterations in function of H
% subplot(1,2,1);
% plot(H,iter, 'red-o', 'LineWidth', 2,MarkerSize= 8)
% title('Evolution of the number of iterations in function of H', 'FontWeight', 'bold')
% xlabel('Length of subdomain H', 'FontWeight', 'bold')
% ylabel('Number of iterations', 'FontWeight', 'bold')
% 
% % computational time in function of the number of subdomains
% subplot(1,2,2);
% plot(rNtab,time, 'red-o', 'LineWidth', 2,MarkerSize= 8)
% title('Evolution of the computational time in function of the number of subdomains', 'FontWeight', 'bold')
% xlabel('Number of subdomain Ns', 'FontWeight', 'bold')
% ylabel('Computational time (s)', 'FontWeight', 'bold')

%% (Q2.6) : Continioning of the assembled preconditioned SP

eigenvals_preSP_SP = abs(eig(preSP_SP));
cond_preSP_SP = max(eigenvals_preSP_SP)/min(eigenvals_preSP_SP);

%% (Q2.7) : Direct resolution of the assembled dual interface problem
tic;
Ub_dual = dualSchur(SDsub, BPsub, Rbsub, blocke, blockAD, Ns, Sd, Bd, G);
elapsed_time_dual = toc;
fprintf('Elapsed time of the dual algorithm: %f seconds\n', elapsed_time_dual);

%% (Q2.8) : Resolution of the assembled dual interface problem using the distributed projected conjugate gradient
% Distributed Projected CG
tic;
[Ub_FETI, k_FETI_dist] = FETI_PCPG(G, SDsub, SPsub, BPsub, Rbsub, ADsub, blockBd, blocke, blockAD, blockSP, Ns, m);
elapsed_time_FETI = toc;
fprintf('Elapsed time of the distributed FETI algorithm: %f seconds\n', elapsed_time_FETI);

% Assembled Projected CG
tic;
[Ub_FETI_bis, k_FETI_assembled] = FETI_PCPG_bis(SDsub, BPsub, Rbsub, G, Sd, Bd, blocke, blockSP ,blockBd, blockAD, Ns, m);
elapsed_time_FETI_assembled = toc;
fprintf('Elapsed time of the assembled FETI algorithm: %f seconds\n', elapsed_time_FETI_assembled);

%% Evaluation of scalability and computational time of FETI
% m = 100;
% iter_dist = []; % Number of iterations
% iter_assembled = [];
% time_dist = []; % Computational time
% time_assembled = [];
% N = 3; % Number of elements (fixed)
% Ntab1= 2:2:200; % Number of subdomains
% H = [];
% rNtab = [];
% for Ns = Ntab1
%     [mapGlobalNodes, interfaceNode, APsub, ADsub, Ksub, SPsub, BPsub, Fsub, SDsub] = decompose_primal_dual(E,S,L,Fd,Ns,N);
%     [Rbsub, blockRb] = rigid_body_mode(SPsub, Ns); % rigid body modes
%     [SP, G, BP, blockAP, blockSP, blockSD, blockBP, blocke, blockBd,blockAD, Sd, Bd] = assemble(SPsub, SDsub, BPsub, Rbsub, APsub, ADsub, Ns);
%     Ub_iter = zeros(Ns-1,1);
%     % Distributed
%     tic;
%     [Ub_FETI, k_dist] = FETI_PCPG(G, SDsub, SPsub, BPsub, Rbsub, ADsub, blockBd, blocke, blockAD, blockSP, Ns, m);
%     t_dist = toc;
%     iter_dist(end+1) = k_dist;
%     time_dist(end+1) = t_dist;
%     % Assembled
%     tic;
%     [Ub_FETI, k_assembled] = FETI_PCPG(G, SDsub, SPsub, BPsub, Rbsub, ADsub, blockBd, blocke, blockAD, blockSP, Ns, m);
%     t_assembled= toc;
%     iter_assembled(end+1) = k_assembled;
%     time_assembled(end+1) = t_assembled;
%     rNtab(end+1) = Ns;
%     H(end+1) = L / Ns; % length of subdomains
% end
% % Create figure figure
% figure
% 
% % number of iterations in function of H
% subplot(2,2,1);
% plot(H,iter_dist, 'red-o', 'LineWidth', 2,MarkerSize= 8)
% title('Number of iterations = f(H), FETI distributed', 'FontWeight', 'bold')
% xlabel('Length of subdomain H', 'FontWeight', 'bold')
% ylabel('Number of iterations', 'FontWeight', 'bold')
% 
% subplot(2,2,2);
% plot(H,iter_assembled, 'red-o', 'LineWidth', 2,MarkerSize= 8)
% title('Number of iterations = f(H), FETI assembled', 'FontWeight', 'bold')
% xlabel('Length of subdomain H', 'FontWeight', 'bold')
% ylabel('Number of iterations', 'FontWeight', 'bold')
% 
% % computational time in function of the number of subdomains
% subplot(2,2,3);
% plot(rNtab,time_dist, 'green-o', 'LineWidth', 2,MarkerSize= 8)
% title('computational time = f(Ns), FETI distributed', 'FontWeight', 'bold')
% xlabel('Number of subdomain Ns', 'FontWeight', 'bold')
% ylabel('Computational time (s)', 'FontWeight', 'bold')
% 
% subplot(2,2,4);
% plot(rNtab,time_assembled, 'green-o', 'LineWidth', 2,MarkerSize= 8)
% title('Computational time = f(Ns), FETI assembled', 'FontWeight', 'bold')
% xlabel('Number of subdomain Ns', 'FontWeight', 'bold')
% ylabel('Computational time (s)', 'FontWeight', 'bold')

