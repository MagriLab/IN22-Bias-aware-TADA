addpath(genpath(pwd)); plot_settings(20);fprintf('\n')
% ======================================================================= %
model = 'Galerkin_low-sqrt';
range_beta  =   1E6;

range_tau   =   1E-4:5E-4:0.1;

fn_bifurcation_diagram(model, range_beta, range_tau)