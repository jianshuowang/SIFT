function [Delta_t, t_error, t, IFT] = analyseSpectrum(lambda,spec,method,debug)
% 
%   Summary:    Function to extract the time delay Delta_t between two laser pulses 
%               by analysing the spectral interference pattern they generate. 
%
%   Inputs:     lambda  - the wavelength array corresponding the the CCD. 
%               spec    - the spectrum of the pulses. 1D spectrum!
%               method  - The method to be used to analyse the spectrum: options: 'fft'
%               debug   - This flag will help to debug any issues
%
%   Outputs:    Delta_t - the time delay between the pulses
%               t_error - Error in the timing measurement. For the 'fft' method, this is
%                         the FWHM of the timing sideband.
%               t       - the time axis from the spectrum's Fourier transform
%               FT      - The Fourier transform of the spectrum
%
%---------------------------------------------------------------
% Author:         Rob Shalloo
% Affiliation:    Imperial College London 
%                 & John Adams Institute for Accelerator Science
% email:          rob.shalloo@gmail.com
% Website:        https://github.com/rob-shalloo
% Created:        November 2018
%---------------------------------------------------------------



% Constants
c = 2.9972e8;        % The speed of light in vacuum

% The data is to be Foruier transformed and thus will need a time axis to
% extract the delay between the pulses. First we convert the horizontal
% axis from wavelength to angular frequency
omega = 2*pi*c./lambda;

% To convert this to a time axis, we need to figure out the sampling
% frequency. YOU NEED TO BE VERY CAREFUL WHEN DEFINING THE SAMPLING FREQUENCY - this
% is because it changes with the frequency due to the nonuniformity of the
% omega axes. Rather than just pick the center value of the sampling
% frequency we will resample the data by interpolating the spectrum onto a
% a new uniform angular frequency axis

% Resample the data onto a uniformly spaced axis
omega_uni = linspace(min(omega),max(omega),length(omega));
spec_uni = interp1(omega,spec,omega_uni);

% Find the new time axis by finding the (now uniform) sampling frequency
tSampleFreq = length(omega_uni)*abs(omega_uni(10)-omega_uni(9));
dt = 2*pi/tSampleFreq;
t = [-dt*length(omega_uni)/2 : dt : dt*length(omega_uni)/2 - dt ];


if method == 'fft'
   % Next step is to take an inverse Fourier transform of the resampled spectrum
   IFT = fftshift(ifft(ifftshift(spec_uni)));
 

   % Find the sideband peak, this corresponds to a first estimate of the 
   % temporal separation between the pulses
   [~,locs] = findpeaks(abs(IFT),t,'SortStr','descend');

   % The timing peak will be the second largest peak (the largest peak in t>0)
   % You could alternatively use the negative peak, but you'd have to be clinically insane
   % To find the peak we simply loop through the peaks in descending order and look for
   % the first positive peak that isn't the largest peak in the time spectrum.
   cntr = 2;
   while locs(cntr) < 0
       cntr = cntr+1;
   end

   Delta_t_estimate = locs(cntr);

   % Now that we have our estimate, reduce our spectrum to only positive values
   [~,indx] = find(t>=0);
   t = t(indx:end);
   IFT = IFT(indx:end);

   % Now hone in on the peak using a fitting method, where we fit two Gaussians to
   % the positive side of the spectrum.
   fitData = fitTimingSpectrum(t, IFT,Delta_t_estimate*1e12,debug);

   % Now extract our measurements from the fit coefficients
   FHWM_centralPeak = fitData(2)*sqrt(2*log(2));
   % This gives an estimate of the measurement error
   t_error  = fitData(3)*sqrt(2*log(2)); % FWHM_timing Peak.
   Delta_t = fitData(4); 
else 
   disp("You've selected a method that Rob hasn't written yet...")
   Delta_t = 0;
   t_error = 0;
end






if debug
   figure
   IFT_debug = abs(IFT);
   plot(t*1e12,abs(IFT_debug),'Color',[.4 .4 1],'LineWidth',2)
   [~,indx] = min(abs(t-Delta_t));
   hold on
   plot(t(indx)*1e12,IFT_debug(indx),'x')
   hold off
   grid on
   xlim([0,2*Delta_t*1e12])
   pbaspect([4 1 1])
end

end



