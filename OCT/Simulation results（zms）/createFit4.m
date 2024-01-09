function [fitresult, gof] = createFit4(Xdata2, Ydata2)
%CREATEFIT4(XDATA2,YDATA2)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : Xdata2
%      Y Output: Ydata2
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  另请参阅 FIT, CFIT, SFIT.

%  由 MATLAB 于 28-Feb-2023 09:50:27 自动生成


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( Xdata2, Ydata2 );

% Set up fittype and options.
ft = fittype( 'gauss1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf 0];
opts.StartPoint = [0.745631568310911 21 1.1537533901855];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData );
% legend( h, 'Ydata2 vs. Xdata2', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% % Label axes
% xlabel( 'Xdata2', 'Interpreter', 'none' );
% ylabel( 'Ydata2', 'Interpreter', 'none' );
% grid on


