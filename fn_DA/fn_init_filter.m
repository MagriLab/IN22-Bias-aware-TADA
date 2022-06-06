function Filter = fn_init_filter(Truth, model_f, model_t)
% Funtion that (1) initialises standard deviations, ensemble
% sizes, etc.; (2) modifies the defined filter time values (initial filter
% time, stop filter time, and time between analysis) to be consistent with
% the simulation times; and (3) initialises the ensemble quantities.
% ----------------------------------------------------------------------- %
    Filter          =   Truth.P;
    Filter.law      =   Truth.P.law;
    % ======================== Filter parameters ======================= %
    Filter.m        =   50;     % Ensemble size
    Filter.sig      =   0.01;    % Std of state initial condition 
    Filter.sig_PE   =   0.15;    % Std of parameters IC
    Filter.sig_mic  =   0.01;   % Std of microphone observations 
    Filter.k_meas   =   40;
    % ================== Type of simulation parameters ================== %

    Filter.E_Params     =   2;
    Filter.E_State      =   1;  
    if Filter.E_Params
        Filter.beta     =   Filter.beta * 1.2;
        if Filter.E_Params == 2
            Filter.tau  =   Filter.tau * .8;
        end
    end
    Filter.E_Bias       =   1;   
    Filter.model_f      =   model_f;
    Filter.model_t      =   model_t;
    Filter.inflation    =   1.01; 
    Filter.N            =   Filter.N + Filter.E_Params;
    % ====================== Filter time parameters ===================== %    
    t_start         =   1;       % Start filtering (after transient) [s]
    t_stop          =   1.5;     % Stop filtering [s]
    t_max           =   2.5;       % Stop simulation [s]
    % Approximate parameters to match the data timestamp
    Filter.t_start	=   approx(Truth.t_mic, t_start);
    Filter.t_stop	=   approx(Truth.t_mic, t_stop);
    Filter.t_max	=   approx(Truth.t_mic, t_max);
    Filter.dt_mic   =   Truth.t_mic(2) - Truth.t_mic(1);
    % ==================== Mics sin/cosine matrices ===================== %    
    c_0     =   Filter.Mean.c_0;
    Filter.sin_omj_mic  =   sin(Filter.omega_j./c_0 * Filter.x_mic)';
    Filter.cos_omj_mic  =   cos(Filter.omega_j./c_0 * Filter.x_mic)';
    % ================== Observation operation matrix =================== %       
    Filter.y_0	=   horzcat(zeros(1,Filter.N), ones(1,Filter.N_mic));
    q       	=   sum(Filter.y_0~=0);    
    M     	=   zeros(q,length(Filter.y_0));   % 
    iq    	=   1;
    for i = 1:length(Filter.y_0)
        if Filter.y_0(i) ~= 0  
            M(iq,i) =   1;      
            iq      =   iq + 1;    
        end
    end
    Filter.M    =   M;
    % =================================================================== %    
    Filter = orderfields(Filter);
end

%% ====================================================================== %
function approx_val = approx(t, val)
    [~,I]       =   min(abs(t-val));
    approx_val  =   t(I);
    
    if abs(approx_val - val) > 1E-3
        error('approximaiton too big, redefine time values')
    end
end