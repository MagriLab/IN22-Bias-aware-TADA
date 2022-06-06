function fn_bifurcation_diagram(model, range_beta, range_tau)
    pks_vec     =   cell(length(range_beta), length(range_tau));
    pks_n_vec   =   cell(length(range_beta), length(range_tau));
    % Setup the simulation
    for betai = 1:length(range_beta)
        Sim                     =   setup_sim(model);    
        Sim.Measurement.Time    =   3;     
        Sim.Mean.beta           =   range_beta(betai);
        for taui = 1:length(range_tau)
            Sim.Mean.Tau    =   range_tau(taui);
            % =========================================================== %    
            if contains(model, 'wave','IgnoreCase',true); fprintf('Wave approach:')      
                % Time march the simulation. This is an origunal code.
                Geom   	=   Sim.Geom;
                Mean	=   Sim.Mean;
                t_end   =   Sim.Measurement.Time;
                [ts,Gs,Hs,~] = fn_time_marching(t_end,Sim.dt,Geom,Mean); 
                % Compute pressure at microphone location at the measurement frequency
                [~, p_mic, ~]   =   GH_to_pu(Sim,ts,Gs,Hs);
            elseif contains(model, 'Galerkin','IgnoreCase',true)
                P   =   setup_P_dim(Sim);        
%                 fn_print_CFL(P);
                [t_vec, psi_vec]	=   ode45(@(t,y) gov_eqns_dim_tau(t,y,P,P.law), ...
                                          (P.t_min:P.dt:P.t_max), P.IC);        
                % Compute pressure and velocity at the microphones' locations
                [~, p_mic, ~]   =   psi_to_pu(t_vec,psi_vec,P);
            else
                error('Type of solver not defined')
            end
            if size(p_mic,1)>size(p_mic,2)
                p_mic = p_mic';
            end
            p_mic   =   p_mic(4,end-2/P.dt:end); % Fourth mic is at flame location
            % =========================================================== % 
            dn = 1;
            X    =   p_mic(1:dn:end);
            X   =   X - mean(X);
            % Check if the data is a fixed point
            x_zeros	=   find(abs(X) < 1e-5);
            if length(x_zeros) == length(X)
                pks_vec{betai,taui} 	=   0;
                pks_n_vec{betai,taui} 	=   0;
            else
                %  Find peaks in the pressure 
                x_pks                   =   X + abs(min(X));
                [pks,~]                 =   findpeaks(x_pks);
                [pks_n,~]               =   findpeaks(-x_pks);
                pks_vec{betai,taui}     =   pks - abs(min(X));
                pks_n_vec{betai,taui}   =   -pks_n - abs(min(X));
            end
        end
    end
    
    betai = length(range_beta);
    taui = length(range_tau);
    [~, sweep_var]   =   max([betai, taui]);
    [~, fixed_var]   =   min([betai, taui]);
    if betai > 4 && taui > 4
        prompt      =   'Sweep beta[1] or tau[2]? ';
        sweep_var	=   input(prompt);
        if sweep_var == 1; fixed_var = 2;
            prompt      =   'Which value of tau in seconds? ';
            tau_ref     =   input(prompt);
            [~, idx]   	=   min(abs(range_tau - tau_ref));
            tau_ref     =   range_tau(idx); 
            disp(['closest value = ',num2str(tau_ref)]);
        else; fixed_var = 1;
            prompt      =   'Which value of beta in Watts/m2*sqrt(m/s)? ';
            beta_ref     =   input(prompt);
            [~, idx]   	=   min(abs(range_beta - beta_ref));
            beta_ref     =   range_beta(idx); 
            disp(['closest value = ',num2str(beta_ref)]);
        end
    end
    params          =   {range_beta, range_tau};
    names           =   {'$\beta', '$\tau'};
    sweep_param     =   params{sweep_var};
    fixed_param     =   params{fixed_var};
    c = 'b.';
    for fi = 1:length(fixed_param)
        fixed   =   fixed_var(fi); 

        figure; hold on
        title([names{fi}, ' = ',num2str(fixed_param(fixed), '%2.2e'),'$'], 'FontSize', 18)
        xlabel([names{sweep_var}, '$'])
        ylabel('$p''_{\mathrm{f}\,\mathrm{max,min}}$')

        if fixed_var == 1
            PKS     =   pks_vec(fi,:);
            PKS_N   =   pks_n_vec(fi,:);
        else
            PKS     =   pks_vec(:,fi);
            PKS_N   =   pks_n_vec(:,fi);
        end
        % Plot peaks in a bifurcation diagram -----------------------------
        for ii = 1:length(sweep_param)
            plot(sweep_param(ii) * ones(1,length(PKS_N{ii})), PKS_N{ii},c,'markersize', 6)
            plot(sweep_param(ii) * ones(1,length(PKS{ii})),PKS{ii},c, 'markersize', 6)
        end
        xlim([sweep_param(1), sweep_param(end)])
    end
    
end