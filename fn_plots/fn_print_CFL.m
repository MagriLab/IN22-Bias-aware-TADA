function CFL = fn_print_CFL(varargin)

    if nargin == 1
        P       =   varargin{1};
        N_c     =   P.N_c;
        dt      =   P.dt;
        tau_adv =   P.tau_advection;
    else
        N_c     =   varargin{1};
        dt      =   varargin{2};
        tau_adv =   varargin{3};
    end

    [~, gc]     =   Cheb(N_c);
    X           =   (gc + 1) ./2;   % Transform g into [0,1]
    dX          =   X(2:end) - X(1:end-1);

    CLF_max     =   1/tau_adv * dt / min(dX);
    CLF_min     =   1/tau_adv * dt / max(dX);
    CLF_mean    =   1/tau_adv * dt / mean(dX);

    fprintf('\n CLF number: max = %2.2f, min = %2.2f, mean = %2.2f\n', ...
            CLF_max, CLF_min, CLF_mean)
    CFL     =   [CLF_max, CLF_min, CLF_mean];
end