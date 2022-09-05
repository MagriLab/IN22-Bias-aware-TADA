close all; clear all; clc
addpath(genpath(pwd)); plot_settings(20);fprintf('\n')

% % Uncomment this if you get errors like "Unable to resolve py.module"
% terminate(pyenv)
% pyenv('Version', '/usr/bin/python3.8'); % change the second argument to
% the location of your python executable
% pyenv("ExecutionMode","OutOfProcess"); % make sure this is the case

% Select data folder
folder = 'DATA';
if ~isfolder(folder)
    mkdir(folder); 
end

%% ======================== SELECT WORKING MODELS ======================= %
% Defined as a two-part string. First, the desired solver, this can be the
% Galerkin method ('Galerkin') or the wave approach, using José Aguilar
% codes, with either high ('Wave_high') low ('Wave_low') order setting. The
% second par defines the heat release law. This can be the square-root
% model ('sqrt'), the modified square-root model ('sqrt_mod'), the 
% arc-tangent model ('tan'), or the G-equation('G_eqn').
% ----------------------------------------------------------------------- %
model_true      =   'Wave_high-G_eqn';
model_forecast	=   'Galerkin_low-sqrt';

% -------------- Select IC for parameters in forecast model ------------- %
train_beta      =   5E6;
train_tau       =   1.5E-3;

% --------------------------- Compare models ---------------------------- %
% Higher fidelity model (truth)
name  =   [folder,'/Truth_', model_true];
try load(name); fprintf(['Loading true data: ',model_true,'\n'])
catch; fprintf(['Creating true data: ',model_true,'\n'])
    Sim_t   =   fn_create_observations(model_true); 
    save([name,'.mat'], 'Sim_t'); % plot_xf(Sim_t);
end
% Forecast model comparison with same parameters
model_forecast  =   [model_forecast,'_', num2str(train_beta,'%1.1e'),...
                    '_', num2str(train_tau,'%1.1e')];
name            =   [folder,'/Truth_', model_forecast];
try load(name); fprintf(['Loading forecast data: ',model_forecast,'\n'])
catch; fprintf(['Creating forecast solution: ',model_forecast,'\n'])
    Sim_f = fn_create_observations(model_forecast, train_beta, train_tau); 
    save([name,'.mat'], 'Sim_f'); % plot_xf(Sim_f) % Plot models
end 
% ------------------------ Compute and plot bias ------------------------ %
% specify the filename as third argument to save the data
models      =	[model_true,'_+_', model_forecast];
[B, t_B]    =   fn_plot_bias(Sim_t, Sim_f, [folder,'/BIAS_',models]); 

%% =========================== LOAD ESN DATA ============================ %
name   =   [folder,'\ESN_',models];
try ESN_P = load(name); fprintf('Loading ESN data \n')
catch; fprintf('Creating ESN data...\n')
    pyrunfile('ESN_RVC_validation.py', filename=models, folder=folder,...
            BETA = train_beta, TAU = train_tau); fprintf('done!\n')      
    ESN_P   =   load(name);
end

% ----------------- Check ESN comparing to actual bias ------------------ %
fn_check_esn(Sim_t, Sim_f, B, t_B, ESN_P); pause(.5);
    

%% ============================ ASSIMILATION ============================ %
rng(2)  % For reproducibility
% ------------------------- Initialise filter --------------------------- %
Filter          =   fn_init_filter(Sim_f, model_forecast, model_true);

% ------------------- Define the time of observations ------------------- %
[obs_idx, Filter] = fn_define_observation_indices(Filter, ESN_P);

% ----------------------- Apply Data Assimilation ----------------------- %
Truth.t_mic     =   Sim_t.t_mic(obs_idx);
Truth.p_mic     =   Sim_t.p_mic(:,obs_idx)';

Filter.model_t  =   model_true;
print_assimilation_type(Filter, ESN_P)

if contains(model_forecast, 'Galerkin')
    [t_KF,Aa_KF,Obs,t_U,U,to,Uo] = EnSRKF_Galerkin(Truth, Filter, ESN_P);
else 
    error('EnSRKF_Wave not yet developed. Select Galerkin forecast model')
end


%% ============================ PLOT RESULTS ============================ %
fn_plot_results(Filter, ESN_P, Sim_t, Truth, t_KF, Aa_KF, ...
                Obs, t_U, U, to, Uo, models)

