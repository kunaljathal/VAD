% *************************************************************************
% Kunal Jathal
% N19194426
% DST 2 - Final Project
% 
% Name:     Periodicity
%
% Purpose:  Helper Function
%
% Description:
%
% This function computes the periodicity of a given signal. It does so by
% making use of a special modified version of the normalized
% autocorrelation function:
%
% P(h) = (ACF(h)/(n - h)) / (ACF(0)/n)              Tmin <= h <= Tmax
% 
% where ACF is the autocorrelation function given by
%
% ACF(h) = (1/n) Sum(t=1 to n - h) X(t+h)*X(t)      Tmin <= h <= Tmax
%
% Essentially, h is the lag value that varies from Tmin to Tmax, which are
% the lower and upper limits of the pitch period respectively. We evaluate
% only pitch periods that lie within the range of human voice pitches, and
% the value of h that gives us the maximum P(h) is the pitch period. The
% maximum P(h) value itself is the periodicity of the frame.
% 
% 
% Usage
% 
% peakPeriodicity = Periodicity(input, fs) returns the periodicity of the
% frame 'input' that is sampled at sampling frequency 'fs'.
% 
% *************************************************************************

function peakPeriodicity = Periodicity(input, fs)

% First, we will Center Clip the signal for improved computational
% efficiency. Center clipping involves implementing the following:
% 
% y(n) = clc[x(n)] = x(n) - Cl,   x(n) >= Cl
%                    0,          |x(n)| < Cl
%                    x(n) + Cl,   x(n) <= -Cl
% 
% where Cl = Clip Threshold

duration = length(input);
clipThresholdFactor = 0.5;
clipThreshold = clipThresholdFactor * max(input);
output = zeros(1, duration);

% Center Clipping algorithm
for i=1:length(input)
    if (input(i) >= clipThreshold)
        output(i) = input(i) - clipThreshold;
    elseif (abs(input(i)) < clipThreshold)
        output(i) = 0;
    elseif (input(i) < - clipThreshold)
        output(i) = input(i) + clipThreshold;
    end
end

% Output now contains the center-clipped signal. Let's get it's ACF.
acf = xcorr(output);

% Now, we want the ACF only for lag values that fall within the pitch
% period limits. Let's take 2.5ms to 15ms, i.e. with a sampling rate of
% fs, that corresponds to:
lowerLagMS = 2.5;
upperLagMS = 15;
lowerLag = round((lowerLagMS * fs) / 1000);
upperLag = round((upperLagMS * fs) / 1000);

% The ACF is symmetric around 0, meaning it goes from -n to n (lag values).
% We want the ACF values for lags going from the lower to the upper lag
% limits. Because there are 2n-1 ACF values, this maps to n + the lag we
% want.
acfSectionalized = acf(duration+lowerLag:duration+upperLag);

% We normalizing the ACF with n - h and then energy per sample.
hVector = lowerLag:upperLag;
lagArray = duration - hVector;
acfNormalized = acfSectionalized ./ lagArray;

energyNormalized = acf(duration)/duration;

finalacf = acfNormalized / energyNormalized;

peakPeriodicity = max(finalacf);


return
