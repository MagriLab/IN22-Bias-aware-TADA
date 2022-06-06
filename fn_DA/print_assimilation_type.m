function print_assimilation_type(Filter,Sim_f, ESN_P)
formatSpec = ['\n%% ------------------ Type of assimilation and parameters --------------------%%\t\n',...
            '%% \t Truth model: %s \t Forecast model: %s \t %%\t\n',...
            '%% \t State estimation: %i \t\t\t\t Parameter estimation: %i \t\t\t %%\t\n',...
            '%% \t dt_analysis: %d \t\t\t\t\t Ensemble size: %i \t\t\t\t %%\t\n',...
            '%% \t Mics uncertainty: %2.2f \t\t\t Initial param uncertainty: %1.2f \t %%\t\n',...
            '%% \t Training beta: %2.2e \t\t\t Initialised beta: %2.2e \t\t %%\t\n',...
            '%% \t Time delay tau: %2.2e \t\t\t Inflation factor: %2.2f \t\t\t %%\t\n',...
            '%% \t ESN trainig time: %2.2f \t\t\t Washout time: %2.4f \t\t\t\t %%\t\n',...
            '%%----------------------------------------------------------------------------%%\t\n'];

fprintf(formatSpec, ...
        Filter.model_t, ...
        Filter.model_f, ...
        Filter.E_State, ...
        Filter.E_Params, ...
        Filter.k_meas, ...
        Filter.m, ...
        Filter.sig_mic, ...
        Filter.sig_PE, ...
        Sim_f.P.beta, ...
        Filter.beta, ...
        Filter.tau, ...
        Filter.inflation, ...
        ESN_P.training_time,...
        ESN_P.N_washout*ESN_P.dt...
        )
end