function plot_settings(varargin)
%% Printing options
FontName = 'monospace';
if nargin == 1
    FontSize = varargin{1};
else
    FontSize = 20;
end
AxLineWidth = 0.75;
LineWidth = 1.00;
MarkerSize = 5;

% Setting LaTeX as text interpreter
set(groot,'DefaultTextInterpreter','latex', ...
    'DefaultAxesTickLabelInterpreter','latex', ...
    'DefaultLegendInterpreter','latex');

% Setting several default options for the graphics
set(groot, 'DefaultFigureColor','White' , ...
'DefaultFigurePaperType', 'a4letter', ...
'DefaultAxesColor', 'white', ...
'DefaultAxesFontUnits', 'points',...
'DefaultAxesFontSize', FontSize, ...
'DefaultAxesTitleFontSize', 1.5, ...
'DefaultAxesFontAngle', 'normal', ...
'DefaultAxesGridLineStyle', '-', ...
'DefaultAxesGridAlpha', 0.25, ...
'DefaultAxesGridColor', [0, 0, 0], ...
'DefaultAxesInterruptible', 'on', ...
'DefaultAxesLayer', 'Bottom', ...
'DefaultAxesNextPlot', 'replace', ...
'DefaultAxesUnits', 'normalized', ...
'DefaultAxesXcolor', [0, 0, 0], ...
'DefaultAxesYcolor', [0, 0, 0], ...
'DefaultAxesZcolor', [0, 0, 0], ...
'DefaultAxesBox','on',...
'DefaultAxesVisible', 'on', ...
'DefaultAxesLineWidth', AxLineWidth, ...
'DefaultLineLineWidth', LineWidth, ...
'DefaultLineMarkerSize', MarkerSize, ...
'DefaultTextColor', [0, 0, 0], ...
'DefaultTextFontUnits', 'Points', ...
'DefaultTextFontName', FontName, ...    % Not applied when LaTeX is used
'DefaultAxesFontName', FontName, ...    % Not applied when LaTeX is used
'DefaultTextFontSize', FontSize, ...
'DefaultTextVerticalAlignment', 'middle', ...
'DefaultTextHorizontalAlignment', 'center')


% Setting a position for the figure
set(groot,'DefaultFigurePosition', [360   198   560   420]);

end
