%% Create figure for output
%close all
clear all
fig1 = figure;
c = 3e8;

%% Define the pulses and their separation
tau     = 40e-15;           % pulse duration in femto second FWHMM pulse dur
taue2   = tau/(sqrt(2*log(2)));   % 1/e^2 pulse duration
amp1    = 1;                % Amplitude of pulse one
amp2    = 1;                % Amplitude of pulse two
dt      = 4.5*1000e-15;          % Separation between the pulses
Nt0     = 2^20;             % Number of points in initial time grid (make big)
lambda = 800e-9;
omega = 2*pi*c/lambda;

    % Define the pulse
    t0 = linspace(-20e-12,20e-12,Nt0);
    p1 = sqrt(amp1*exp( - 2*(t0-dt/2).^2/(taue2)^2)).*exp(1i*omega*t0);
    p2 = sqrt(amp2*exp( - 2*(t0+dt/2).^2/(taue2)^2)).*exp(1i*omega*t0).*exp(1i*1.5);
    y0 = p1 + p2;
    %y0 = amp1*exp( - 2*(t0-dt/2).^2/(taue2)^2) + amp2*exp( - 2*(t0+dt/2).^2/(taue2)^2);
    %y0 = sqrt(y0).*exp(1i*omega*t0);
    
    
    % Plot it 
    fig1
    subplot(4,1,1)
    plot(t0,abs(y0).^2,'LineWidth',2,'Color',[0 0 0])
    hold on
    plot(t0,real(y0),'LineWidth',.5,'Color',[.3 .3 .3])
    hold off
    xlim([-dt,dt])

%% Fourier Transform Data and Interpolate onto the measurement grid - your spectrometer
dLambda     = .02e-9;       % Spectral resolution of the spectrometer
N           = 4096;         % Number of pixels in the spectrometer (horizontal)
lambda0     = 800e-9;       % Center wavelength of spectrometer

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

%%%%%%%%
%IFT = fftshift(ifft(ifftshift(abs(FT).^2)));
%subplot(4,1,3)
%plot(t0,abs(IFT),'LineWidth',2,'Color',[0 0 0]);

%%%%%%%%


    % Plot it
    subplot(4,1,2)
    hold off
    plot(xFreq,abs(FT).^2,'LineWidth',2,'Color',[0 0 0]);
    xlim([omega-2e14,omega+2e14])
    
    
    % Interpolate onto grid
    lambda = linspace(lambda0- (N/2)*dLambda ,lambda0 + (N/2-1)*dLambda,N);
    spec = interp1(2*pi*c./(xFreq),abs(FT).^2,lambda);
    subplot(4,1,3)
    plot(lambda*1e9,spec,'LineWidth',2,'Color',[0 0 0]);
    xlim([760,840])
    
% Make sure this measured spectrum matches out original.
subplot(4,1,2)
hold on
omeg_meas = 2*pi*c./lambda;
plot(omeg_meas,spec,'.-','LineWidth',.5,'Color',[1,.3,.3])
    
%% Now Fourier transform back

 [Delta_t,FWHM_timingPeak,t,IFT] = analyseSpectrum(lambda,spec,1);
 figure(fig1)
 subplot(4,1,4)
 hold off
 plot(t*1e12,abs(IFT)/max(abs(IFT)),'LineWidth',2,'Color',[0 0 0]);
 xlim([0,2*dt*1e12])
 
 %% Checks and balances
delta_t =  fwhm(t0,abs(y0).^2);
delta_lambda =  fwhm(lambda,spec);
disp(sprintf('Time BandWidth Product = %f',delta_t*delta_lambda/lambda0/lambda0*c));

disp(sprintf('The measured temporal Delay is (%f +/- %f) ps',Delta_t*1e12,FWHM_timingPeak*1e12))