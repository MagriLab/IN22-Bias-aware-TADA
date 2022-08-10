plot_settings(16)
model_forecast	=   'Galerkin_low-sqrt'; 
Sim             =   setup_sim(model_forecast);
Sim.Mean.Tau    =   0.002;

P   =   fn_extract_Sim_properties(Sim, []);    

for Nc  =   [120]
    P.N_c       =   Nc;                     % Num Chebishev modes
    P           =   properties(P);
    %%
    
    CFL = fn_print_CFL(P);
    
    [t, psi]	=   ode45(@(t,y) gov_eqns_dim_tau(t,y,P,P.law), ...
                                      (P.t_min:P.dt:P.t_max), P.IC);
    
    %
    N_m         =   P.N_m;
    N_c         =   P.N_c;
    dt          =   P.dt;
    [Dc, gc]    =   Cheb(N_c);
    X           =   (gc(2:end) + 1) ./2;   % Transform g into [0,1]
    
    %
    advection_modes     =   psi(:,2*N_m+1:2*N_m+N_c);
    u_xf                =   P.cos_omjxf * psi(:,1:N_m)';
    
    
    %%
    figure("Units","normalized",'Position',[0 .2 1 0.5]); tiledlayout(2,5); 
    tl = 0;
    for delay = [0, Sim.Mean.Tau, P.tau_advection]; tl = tl + 1; % 
        
        delay_dt    =   delay / dt;             % [num of dt]
        
        
        X_delay                 =   delay/P.tau_advection;
        u_delay_advection       =   interp1(X, advection_modes', X_delay,"pchip","extrap");
        u_delay                 =   u_xf * 0;
        u_delay(delay_dt+1:end) =   u_xf(1:end-delay_dt);
        
        
        % figure; plot(u_xf)
        %
        nexttile(tl); hold on
        if tl == 2
        title(['$N_c=$',num2str(P.N_c),', CFL$_\mathrm{max} = $ ', num2str(CFL(1)),', CFL$_\mathrm{min} = $ ', num2str(CFL(2)), ', CFL$_\mathrm{mean} = $ ', num2str(CFL(3))])
        end
        plot(t, u_delay_advection, 'DisplayName', 'advection sol.')
        plot(t, u_delay, 'DisplayName', 'actual sol.')
        ylabel(['$u''(t-$',num2str(delay/Sim.Mean.Tau),'$\tau$)'])
        legend()
        xlim([t(end-500), t(end)])
        nexttile(tl+5); hold on
        title(['$\Delta X_\mathrm{interp}= $ ',num2str(min(abs(X-X_delay)),4)])
        plot(t, u_delay_advection- u_delay, 'DisplayName', 'advection sol.')
        xlabel('$t$'); ylabel('$\Delta u''$')
        xlim([t(1), t(end)])
    end
    
    ti = length(t)-2000:length(t);
    nexttile([2,2])
    waterfall(X, t(ti), advection_modes(ti,:)); colormap('parula')
    xlabel('X'); ylabel('$t$ [s]'); zlabel('$v$')
    title(['$N_c=$',num2str(P.N_c),', $\tau=$', num2str(Sim.Mean.Tau), ...
            ' s', ', $\tau_\mathrm{advection}=$', num2str(P.tau_advection), ' s' ])
    hold on
    x_ti1 = 0:dt/P.tau_advection:1;
    v_ti1  =   u_xf(ti(1):-1:ti(1)-length(x_ti1)+1);
    plot3(x_ti1, x_ti1*0 + t(ti(1)), v_ti1, 'k-.','LineWidth',3)
    
    x_ti1 = 0:dt/P.tau_advection:1;
    v_ti1  =   u_xf(ti(end):-1:ti(end)-length(x_ti1)+1);
    plot3(x_ti1, x_ti1*0 + t(ti(end)), v_ti1, 'k-.','LineWidth',3)
    ylim([t(ti(1)) t(ti(end))])
end
    
    
%%  -------------------------------------------------------------------- %%
function P=properties(P)
    % ======================== STATE VECTOR SIZE ======================== %
    P.N_m       =   10;                     % Num acoustic modes
    % P.N_c       =   10;                     % Num Chebishev modes
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
    P.tau_advection     =   0.100;
    if P.tau_advection < P.tau
        P.tau_advection = P.tau;
        disp('advection tau < acustic tau. Forced advection tau = acustic tau')
    end

end