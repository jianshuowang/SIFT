%   Summary:    Test script to investigate how the extracted temporal
%               resolution depends upon the spectral resolution of the 
%               chosen spectrometer        
%       
%---------------------------------------------------------------
% Author:         Rob Shalloo
% Affiliation:    Imperial College London 
%                 & John Adams Institute for Accelerator Science
% email:          rob.shalloo@gmail.com
% Website:        https://github.com/rob-shalloo
% Created:        November 2018
%---------------------------------------------------------------

addpath('../')

%% Create figure for output
close all
clear all

% Constants
c = 2.9972e8;                   % The speed of light in vacuum

% Define the pulses and their separation
tau     = 40e-15;               % pulse duration in femto second FWHMM pulse dur
taue2   = tau/(sqrt(2*log(2))); % 1/e^2 pulse duration
amp1    = 1;                    % Amplitude of pulse one
amp2    = 1;                    % Amplitude of pulse two
Nt0     = 2^20;                 % Number of points in initial time grid (make big)
dt      = 5e-12;               % Separation between the pulses
lambda  = 800e-9;               % fundamental wavelength of the pulse
omega   = 2*pi*c/lambda;        % Convert to angular frequency


% Define Spectrometer Parameters
dLambdas     = linspace(.0001e-9,.4e-9,200); % Spectral resolution of the spectrometer
lambda0     = 800e-9;                        % Center wavelength of spectrometer
requiredBandwidth = 60e-9;                   % Define the bandwidth required on the spectrometer
                                             % we will calculate the required number of pixels to achieve
                                             % this on every loop and use that value. 

% Create empty arrays for the data to be calculated. 

t_meas = zeros(1,length(dLambdas));             
t_meas_err = zeros(1,length(dLambdas));

% Loop through the differing values of spectral resolution. 
for i = 1:length(dLambdas)
    
    % we want to ensure that we cover a minimum bandwidth so calculate the necessary num of pixels
    dLambda = dLambdas(i);
    N = ceil(requiredBandwidth/dLambda);  

    % Define the pulse
    t0 = linspace(-20e-12,20e-12,Nt0);
    p1 = sqrt(amp1*exp( - 2*(t0-dt/2).^2/(taue2)^2)).*exp(1i*omega*t0);
    p2 = sqrt(amp2*exp( - 2*(t0+dt/2).^2/(taue2)^2)).*exp(1i*omega*t0).*exp(1i*1.5);
    y0 = p1 + p2;

    % Fourier Transform the temporal profile of the pulses
    FT = fftshift(fft(ifftshift(y0)));
    xSampleFreq = length(t0)/(t0(end)-t0(1));
    xFreq = [-xSampleFreq/2 : xSampleFreq/Nt0 : (xSampleFreq/2-xSampleFreq/Nt0) ];
    xFreq = 2*pi*xFreq;

    % Now restrict outselves to postive values only
    [~,indxs] = find(xFreq > 0);
    % and add a small delta to ensure we're always over the line (i.e. no dividing by zero)
    del = 100;
    xFreq = xFreq(min(indxs)+del:end);
    FT = FT(min(indxs)+del:end);
    
    % Interpolate onto the pixels of the chosen spectrometer 
    lambda = linspace(lambda0- (N/2)*dLambda ,lambda0 + (N/2-1)*dLambda,N);
    spec = interp1(2*pi*c./(xFreq),abs(FT).^2,lambda);

    % Now make a measurement of the pulse timing
    [Delta_t,t_error,t,IFT] = analyseSpectrum(lambda,spec,'fft',0);
    t_meas(i) = Delta_t;
    t_meas_err(i) = t_error;
end

% Plot the results!
fig1 = figure
errorbar(dLambdas*1e9,t_meas*1e12,t_meas_err*1e12)
hold on
grid on
plot(dLambdas*1e9,dt*1e12*ones(length(dLambdas)))
xlabel('Spectral Resolution (nm)')
ylabel('"Measured" Pulse Separation (ps)')

fig2 = figure
errorbar(dLambdas/lambda0,t_meas*1e12,t_meas_err*1e12)
hold on
grid on
plot(dLambdas/lambda0,dt*1e12*ones(length(dLambdas)))
xlabel('Spectral Resolution (norm.)')
ylabel('"Measured" Pulse Separation (ps)')
