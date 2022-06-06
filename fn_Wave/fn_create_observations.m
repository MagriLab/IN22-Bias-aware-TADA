function T = fn_create_observations(varargin)
    model = varargin{1};
    if contains(model, 'wave','IgnoreCase',true); fprintf('Wave approach:')  
        % =============================================================== %        
        % Setup the simulation
        Sim     =   setup_sim(model);
        % Time march the simulation. This is an origunal code.
        Geom   	=   Sim.Geom;
        Mean	=   Sim.Mean;
        t_end   =   Sim.Measurement.Time;
        [ts,Gs,Hs,Xis,Qs] = fn_time_marching(t_end,Sim.dt,Geom,Mean); 
        % Compute pressure at microphone location at the measurement frequency
        [t_mic, p_mic, u_mic] = GH_to_pu(Sim,ts,Gs,Hs);

        % Store true data
        T.p_mic         =   p_mic;
        T.u_mic         =   u_mic;
        T.t_mic         =   t_mic;
        T.Gs            =   Gs;
        T.Hs            =   Hs;
        T.ts            =   ts;
        T.Xis           =   Xis;
        T.Qs            =   Qs;
        T.Geom          =   Geom;
        T.Mean          =   Mean;
        T.Measurement	=   Sim.Measurement;
        T.Forcing       =   Sim.Forcing;
        % =============================================================== %
    elseif contains(model, 'Galerkin','IgnoreCase',true)
        % =============================================================== %        
        
        beta = varargin{2};
        tau = varargin{3};
        
        Sim =   setup_sim(model, beta, tau);
        P   =   setup_P_dim(Sim);
        
        [t_vec, psi_vec]	=   ode45(@(t,y) gov_eqns_dim(t,y,P,P.law), ...
                                  (P.t_min:P.dt:P.t_max), P.IC);

        
        % Compute pressure and velocity at the microphones' locations
        [t_mic, p_mic, u_mic] = psi_to_pu(t_vec,psi_vec,P);
        
        if size(psi_vec,1) < size(psi_vec,2)
            psi_vec = psi_vec';
        end
        
        % Output
        T.t     =   t_vec';
        T.psi	=   psi_vec;
        T.P     =	orderfields(P);        
        T.p_mic	=   p_mic';
        T.u_mic	=   u_mic';
        T.t_mic =   t_mic;
        % =============================================================== %
    else
        error('Type of solver not defined')
    end
    T	=	orderfields(T);
end

% ======================================================================= %




% ----------------------------------------------------------------------- %
