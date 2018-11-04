function fitData = fitTimingSpectrum(t, IFT,Delta_t_estimate,debug)
% 
%   Summary:    Fit two gaussians to the timing spectrum to extract the temporal
%               separation of the two pulses 
%
%   Inputs:     t        - the time axis from the spectrum's Fourier transform
%               IFT      - The Inverse Fourier transform of the spectral data
%               Delta_t  - the best estimate for the time delay between the pulses
%               debug    - standard debug flag; will plot the fit.
%
%   Outputs:    fitData  - The coefficients of the fit xData;
%                          b  - amplitude of the sideband relative to central peak
%                          c1 - 1/e^2 radius of DC peak
%                          c2 - 1/e^2 radius of sideband peak
%                          t0 - the temporal separation between the pulses             
%               
%
%---------------------------------------------------------------
% Author:         Rob Shalloo
% Affiliation:    Imperial College London 
%                 & John Adams Institute for Accelerator Science
% email:          rob.shalloo@gmail.com
% Website:        https://github.com/rob-shalloo
% Created:        November 2018
%---------------------------------------------------------------


% Rescale the data!
t = t*1e12;
IFT = abs(IFT/max(IFT));


[xData, yData] = prepareCurveData( t, IFT );

% Get the fit parameters ready to go
ft = fittype( 'exp( -2*x.^2/(c1^2)) + b*exp( -2*(x-t0).^2/(c2^2))', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 0 0 Delta_t_estimate/2];
opts.StartPoint = [1 Delta_t_estimate/4 Delta_t_estimate/4 Delta_t_estimate];
opts.Upper = [1 Delta_t_estimate/2 Delta_t_estimate/2 1.5*Delta_t_estimate];

% Fit the data
[fitresult, gof] = fit( xData, yData, ft, opts );


% The fit result needs to be rescaled back to scale of the original data
fitData = coeffvalues(fitresult);
fitData(1) = fitData(1)*max(IFT);
fitData(2) = fitData(2)/1e12;
fitData(3) = fitData(3)/1e12;
fitData(4) = fitData(4)/1e12;

if debug
    % Plot fit with data.
    figure( 'Name', 'untitled fit 1' );
    h = plot( fitresult, xData, yData );
    legend( h, 'IFT vs. t', 'untitled fit 1', 'Location', 'NorthEast' );
    % Label axes
    xlabel t
    ylabel IFT
    grid on
    hold on
    t2 = linspace(0,2*fitData(4)*1e12,1000);
    plot(t2,exp( -2*t2.^2/((fitData(2)*1e12)^2)) + fitData(1)/max(IFT)*exp( -2*(t2-fitData(4)*1e12).^2/((fitData(3)*1e12)^2) ))
end


