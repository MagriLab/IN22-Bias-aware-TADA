function P  =   setup_P_dim(Sim)
    % =================================================================== %
    % Use Jose's code for defining the mean flow properties. Then extract
    % the averaged values and predefined heat release rate law parameters
    P   =   fn_extract_Sim_properties(Sim, []);    
    % ======================== STATE VECTOR SIZE ======================== %
    P.N_m       =   10;                     % Num acoustic modes

    P.tau_advection     =   0.05; %10 * P.tau;
    if P.tau_advection < P.tau
        P.tau_advection = P.tau;
        disp('advection tau < acustic tau. Forced advection tau = acustic tau')
    end
    P.N_c       =   50;                     % Num Chebishev modes
%     fn_print_CFL(P)



    P.N       	=   P.N_m * 2 + P.N_c;     	% State vector size
    % ========================= MODES FREQUENCY ========================= %
    P.omega_j	=   fn_obtain_omega(P)'; 	% Acoustic modes freqs  [Hz]
    [P.D_c, P.gc] =   Cheb(P.N_c);            % Chebyshev matrix 
    P.X         =   (P.gc + 1) / 2;           % Chebyshev pts [0, 1]
    j           =   (1:P.N_m)';
    P.jpi       =   j * pi;
    c_0         =   P.Mean.c_0;
    P.cos_omjxf = 	cos(P.omega_j'./c_0 * P.x_f);	% Cosine matrix at xf
    P.sin_omjxf	=   sin(P.omega_j'./c_0 * P.x_f);	% Sine matrix at xf
    % ============================= DAMPING ============================= %
    P.C1        =	0.1;                    % Damping parameter 1
    P.C2        =   0.06;                	% Damping parameter 2
    P.zeta      =   P.C1*j.^2 + P.C2*sqrt(j);         % Damping
    % ======================== INITIAL CONDITIONS ======================= %
    P.pert  =   0.05;                       % Modes initial condition
    IC_eta  =   P.pert * ones(P.N_m,1);
    IC_v    =   zeros(P.N_c,1);
    P.IC   	=   [IC_eta; IC_eta; IC_v];
    if contains(P.law, 'tan')
        IC_w    =   zeros(P.N_c,1);
        P.IC   	=   [P.IC; IC_w];
    end
    % =================================================================== %

end