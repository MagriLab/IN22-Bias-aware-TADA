function print_assimilation_type(Filter, ESN_P)
formatSpec = ['\n%% ------------------ Type of assimilation and parameters --------------------%%\t\n',...
            '%% \t Truth model: %s \t\t Forecast model: %s \t %%\t\n',...
            '%% \t State estimation: %i \t\t\t\t Parameter estimation: %i \t\t\t %%\t\n',...
            '%% \t dt_analysis: %d \t\t\t\t\t Ensemble size: %3d \t\t\t\t %%\t\n',...
            '%% \t Number of mics: %d \t\t\t\t\t Inflation factor: %2.2f \t\t\t %%\t\n',...
            '%% \t Mics uncertainty: %2.2f \t\t\t Initial param uncertainty: %1.2f \t %%\t\n',...
            '%% \t Training beta: %2.2e \t\t\t Initialised beta: %2.2e \t\t %%\t\n',...
            '%% \t Training tau: %2.2e \t\t\t Initialised tau: %2.2e \t\t\t %%\t\n',...
            '%% \t ESN trainig time [s]: %2.2f \t\t Washout time [s]: %2.4f \t\t\t %%\t\n',...
            '%% \t tau advection: %2.2e \t\t\t Chebyshev modes: %i \t\t\t\t %%\t\n',...
            '%%----------------------------------------------------------------------------%%\t\n'];
model_f = Filter.model_f;
if contains(model_f, 'e+')
    model_f = model_f(1:end-16);
end
fprintf(formatSpec, ...
        Filter.model_t, ...
        model_f, ...
        Filter.E_State, ...
        Filter.E_Params, ...
        Filter.k_meas, ...
        Filter.m, ...
        Filter.N_mic, ...
        Filter.inflation, ...
        Filter.sig_mic, ...
        Filter.sig_PE, ...
        ESN_P.training_beta, ...
        Filter.beta, ...
        ESN_P.training_tau, ...
        Filter.tau, ...
        ESN_P.training_time,...
        ESN_P.N_washout*ESN_P.dt,...
        Filter.tau_advection,... 
        Filter.N_c...
        )
end