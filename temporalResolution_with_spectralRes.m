% Temporal Resolution Tests
%Varying time bewteen pulses

%% Create figure for output
close all
clear all
fig1 = figure;
c = 3e8;

%% Define the pulses and their separation
tau     = 40e-15;           % pulse duration in femto second FWHMM pulse dur
taue2   = tau/(sqrt(2*log(2)));   % 1/e^2 pulse duration
amp1    = 1;                % Amplitude of pulse one
amp2    = 1;                % Amplitude of pulse two
Nt0     = 2^20;             % Number of points in initial time grid (make big)
dt      = 10e-12;            % Separation between the pulses
lambda = 800e-9;
omega = 2*pi*c/lambda;


%% Define Spectrometer Parameters
dLambdas     = linspace(.0001e-9,.4e-9,200); % Spectral resolution of the spectrometer
t_meas = zeros(1,length(dLambdas));
t_meas_err = zeros(1,length(dLambdas));

for i = 1:length(dLambdas)
    
    % we want to ensure that we cover a min 60nm bandwidth
    dLambda = dLambdas(i);
    N = ceil(60e-9/dLambda);  
    lambda0     = 800e-9;       % Center wavelength of spectrometer

    % Define the pulse
    t0 = linspace(-20e-12,20e-12,Nt0);
    p1 = sqrt(amp1*exp( - 2*(t0-dt/2).^2/(taue2)^2)).*exp(1i*omega*t0);
    p2 = sqrt(amp2*exp( - 2*(t0+dt/2).^2/(taue2)^2)).*exp(1i*omega*t0).*exp(1i*1.5);
    y0 = p1 + p2;
    %y0 = amp1*exp( - 2*(t0-dt/2).^2/(taue2)^2) + amp2*exp( - 2*(t0+dt/2).^2/(taue2)^2);
    %y0 = sqrt(y0).*exp(1i*omega*t0);
    
    FT = fftshift(fft(ifftshift(y0)));
    xSampleFreq = length(t0)/(t0(end)-t0(1));
    xFreq = [-xSampleFreq/2 : xSampleFreq/Nt0 : (xSampleFreq/2-xSampleFreq/Nt0) ];
    xFreq = 2*pi*xFreq;

    % Now restrict outselves to postive values only
    [~,indxs] = find(xFreq > 0);
    % and add a small delta to ensure we're always over the line
    del = 100;
    tmpXFreq = xFreq;
    xFreq = xFreq(min(indxs)+del:end);
    %lam_int = 2*pi*c./(xFreq);
    tmpFT = FT;
    FT = FT(min(indxs)+del:end);
    
    % Interpolate onto grid
    lambda = linspace(lambda0- (N/2)*dLambda ,lambda0 + (N/2-1)*dLambda,N);
    spec = interp1(2*pi*c./(xFreq),abs(FT).^2,lambda);
    [Delta_t,FWHM_timingPeak,t,IFT] = analyseSpectrum(lambda,spec,0);
    t_meas(i) = Delta_t;
    t_meas_err(i) = FWHM_timingPeak;
end

figure
errorbar(dLambdas*1e9,t_meas*1e12,t_meas_err*1e12)
hold on
grid on
plot(dLambdas*1e9,dt*1e12*ones(length(dLambdas)))
xlabel('Spectral Resolution (nm)')
ylabel('"Measured" Pulse Separation (ps)')

figure
errorbar(dLambdas/lambda0,t_meas*1e12,t_meas_err*1e12)
hold on
grid on
plot(dLambdas/lambda0,dt*1e12*ones(length(dLambdas)))
xlabel('Spectral Resolution (norm.)')
ylabel('"Measured" Pulse Separation (ps)')