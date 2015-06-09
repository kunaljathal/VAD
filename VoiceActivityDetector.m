% *************************************************************************
% Kunal Jathal
% N19194426
% DST 2 - Final Project
% 
% VOICE ACTIVITY DETECTOR
%
% Name:     VoiceActivityDetector
%
% Purpose:  Main Function
%
% Description:
%
% This function implements a rad VAD. It works by using voiced band energy
% ratios, periodicity measures, running and dynamic minimum and maximum
% RMS energy estimates, adaptive and noise resistant threshold computation,
% and hangover smoothing to fade in and out speech/noise boundaries. The
% final output is a reconstructed voice only audio file, together with a
% visual representation showing the VAD marked sections and the voice plot
% 
% Usage
% 
% VoiceActivityDetector(input, noise) - Performs VAD on the input audio
% signal 'input'. 'noise' is an optional wav of noise to mix into the
% input audio file.
% 
% For input, the following voice audio files are present in the folder:
% 
% male1.wav     -       sample male voice
% male2.wav     -       sample male voice
% male3.wav     -       sample male voice
% female1.wav   -       sample female voice
% female2.wav   -       sample female voice
% 
% For noise, the following files are present in the folder:
% 
% 'white.wav'   -       standard white noise, constant level
% 'pink.wav'    -       standard pink noise, constant level
% 'babble.wav'  -       generic people talking, dynamically changing levels
% 'music.wav'   -       music used as noise, dynamically changing levels
% 
% The noise file is optional. You can simply provide your own already noisy
% signal directly as 'input' if you'd like, or you can also provide your own
% individual voice and noise files separately.
% 
% NOTE: If you are using only a single noisy input file, please make
% sure the noise is reasonably softer than the voice input. The VAD works
% best when the noise is around 10% of the signal amplitude. Also, ensure
% the first ~300ms of the audio is only noise (no voice).
% 
% Finally , to play around with the accuracy of the VAD, you can try tweaking
% the following variables directly in the code itself:
% 
% frameLength           -   ~10 to 20ms (if you want to go below 10, you 
%                                        will need to change the upperLagMS 
%                                        value in Periodicity.m to be lower)
%                                          
% hangOverThreshold     -   ~3 to 10    (# of unvoiced frames to pad  
%                                        unconditionally after the end of a 
%                                        voiced frame
% 
% noiseScalingFactor    -   ~0.05 to ~0.2 (amplification of noise compared
%                                          to voice signal expressed as a 
%                                          percentage
% 
% *************************************************************************

function VoiceActivityDetector(input, noise)

% Get the voice and noise signals ready
voiceSignal = [];
noiseSignal = [];
inputSignal = [];
fs = 0;
fs1 = 0;
fs2 = 0;

% Noise amplification ratio compared to voice expressed as a %
noiseScalingFactor = 0.1;

% Read the input signal
[voiceSignal, fs1] = wavread(input);
voiceSignal = processAudio(voiceSignal);

if (nargin > 1)
    % Read the noise signal
    [noiseSignal, fs2] = wavread(noise);
    
    % Resample if necessary
    if (fs1 > fs2)
        noiseSignal = upSample(noiseSignal, fs2, fs1);
        fs = fs1;
    elseif (fs2 > fs1)
        voiceSignal = upSample(voiceSignal, fs1, fs2);
        fs = fs2;
    else
        fs = fs1;
    end            
    
    % Normalize noise, etc.
    noiseSignal = processAudio(noiseSignal);

    % Trim noise signal to voice signal length if necessary
    if (length(noiseSignal) > length(voiceSignal))
        noiseSignal = noiseSignal(1:length(voiceSignal));
    end
    
    % Scale noise signal down
    noiseSignal = 0.1 * noiseSignal;
    
    % Mix voice and noise
    inputSignal = voiceSignal + noiseSignal;
else
    % Input signal is already noisy
    inputSignal = voiceSignal;
    fs = fs1;
end

inputSignal = processAudio(inputSignal);

% Play the current voice+noise mix
% sound(inputSignal, fs);

% Get the VAD Marker and the output signal ready
vadMarker = [];
outputSignal = [];

% Choose a frame length in ms. This value can be tweaked and results in
% varying performance. I've found that 20ms works well for most signals.
% This value cannot go lower than the 'upperLagMS' variable in
% Periodicity.m, so if you'd like to use a smaller frame size, please
% change upperLagMS accordingly.
frameLength = 20;
frameLengthSamples = ceil((frameLength * fs)/1000);

% Segment the signal into frames
segmentedSignal = frameSplit(inputSignal, frameLengthSamples);
numFrames = size(segmentedSignal, 1);

% Initialize VAD flag arrays
iVad = zeros(numFrames, 1);
fVad = zeros(numFrames, 1);

% Initialize hang over smoothing variables
hangoverThreshold = 5;
inactiveFrameCounter = hangoverThreshold + 1;
firstVoicedFlag = 0;
applyFadeOut = 0;
bufferVoicedFrames = [];
bufferFilling = 0;
bufferLimit = 4 ;
bufferFull = 0;

% Initialize Energy Variables
energyMin = 0;
energyMax = 0;
deltaEmin = 1;
deltaEmax = 1;
lambda = 1;
threshold = 0;
initialEnergyMin = 0;
initialNoiseDuration = 300;
lowerEnergyVoiceBand = 2000;
upperEnergyVoiceBand = 4000;

% Initialize Periodicity and Energy Ratio thresholds
periodicityThreshold = 1;
energyRatioThreshold = 10;

% Initialize RMSE, periodicity, and energy ratio arrays (debugging use)
RMSEArray = zeros(numFrames, 1);
periodicityArray = zeros(numFrames, 1);
ratioArray = zeros(numFrames, 1);
thresholdArray = zeros(numFrames, 1);
thresholdMarker = [];

% Set the fade in and fade out times for hangover smoothing
fadeInTime = frameLengthSamples;
fadeOutTime = frameLengthSamples*hangoverThreshold;

% Use a simple linear envelope for the fades
fadeInEnvelope =  linspace(0, 1, fadeInTime);
fadeOutEnvelope = linspace(1, 0, fadeOutTime);

% Great! We now have the input signal split up into frames and are ready. 
% Let us compute the thresholds and enforce the VAD process per frame.
for frameCounter=1:numFrames
    frame = segmentedSignal(frameCounter, :);

    % Get the Root Mean Square Energy of the frame
    RMSE = sqrt(mean(sum(frame.^2)));
    
    % Add to the array of RMS energies (debugging only)
    RMSEArray(frameCounter) = RMSE;
    
    % NOTE: An important assumption necessarily made is that the first
    % 100ms or so of the signal will be noise. So we use the first few
    % frames to calculate and set Emin and Emax values.
    if ((frameCounter * frameLength) < initialNoiseDuration)
        if (RMSE > energyMax)
            % Set Emax
            energyMax = 5 * RMSE;
        end
        
        if (energyMin == 0)
            % Initialize Emin
            energyMin = RMSE;            
            initialEnergyMin = energyMin;
        else
            if (RMSE > energyMin)
                energyMin = RMSE;
            end
        end
        
        % The first few frames are assumed to be noise, so no VAD
        % processing is necessary.
        iVad(frameCounter) = 0;
        fVad(frameCounter) = 0;
        
    else                
        % We are into the signal now! Commence actual VAD processing.
        % Update running Energy threshold estimates (Emax, Emin)
        if (RMSE > energyMax)
            % We use a small delta to gradually decrease Emax to compensate
            % for anomalous spikes in energy
            energyMax = RMSE;
            deltaEmax = 1;
        else
            deltaEmax = 0.999;
        end

        if (RMSE < energyMin)
            if (RMSE == 0)
                energyMin = initialEnergyMin;
            else            
                energyMin = RMSE;
            end

            deltaEmin = 1;
        else
            % We use a small delta scaling factor to prevent complications
            % arising from energy dips (anomalies). This keeps Emin rising
            % at a gradual, minimal rate.
            deltaEmin = deltaEmin * 1.001;
        end

        % Threshold computations. Lambda is the non-linear dynamic
        % coefficient used to compute the threshold in a way that makes it
        % resistant and independent of variations in background noise.
        lambda = 1 - (energyMin/energyMax);
        threshold = ((1 - lambda) * energyMax ) + (lambda * energyMin);

        % Keep track of threshold values (debugging use)
        thresholdArray(frameCounter) = threshold;
        
        
        % Get the periodicity of the frame
        periodicity = Periodicity(frame, fs);

        % Add to peridocity array (debugging)
        periodicityArray(frameCounter) = periodicity;
        
        % Get the ratio of the frequencies above and below 2 kHz in the voice
        % band (0 - 4 kHz). We will then use the ratio of these energies to
        % make decisions around the voicing of the frame.
        window = (hamming(size(frame, 2)))';
        
        % FFT computation & normalization
        fftLength = 2^nextpow2(length(frame));
        theFFT = fft(frame.*window, fftLength);
        fftLength = length(theFFT);
        fftSq = (theFFT).*conj(theFFT);
        fftSqNorm = fftSq/fftLength;

        % Compute the voiced energy band ratio
        energyBelow = sum(fftSqNorm(1:round((lowerEnergyVoiceBand/fs) * fftLength)));
        energyAbove = sum(fftSqNorm(round((lowerEnergyVoiceBand/fs) * fftLength):round((upperEnergyVoiceBand/fs) * fftLength)));
        energyRatio = energyAbove/energyBelow;

        % Add to ratio array (debugging)
        ratioArray(frameCounter) = energyRatio;
        
        % Now that we have all the features extracted, make a decision about
        % the frame being voiced or not. We decide that a frame is voiced
        % if either a) the RMS Energy is above the threshold, or b) the
        % signal shows periodicity within the human voice range of pitches,
        % or c) there is significantly higher energy in the upper range of
        % the voiced band as compared to the lower range, pivoting at 2 kHz
        if ((RMSE > threshold) || (periodicity > periodicityThreshold))
            iVad(frameCounter) = 1;
        elseif (abs(energyRatio) > energyRatioThreshold)
            iVad(frameCounter) = 1;
        else
            iVad(frameCounter) = 0;
        end
    end
    
    
    % Before constructing the output signal, we perform hang-over smooting
    % to ensure audibly pleasing transitions from speech to noise and vice
    % versa
    if (iVad(frameCounter) == 0)

        % If the buffer is being filled and we encountered an unvoiced
        % frame, it represents an anomalous spike in energy, so record that
        % and set the limit to as to wipe the buffer out later and put in
        % an unvoiced frame later instead of a voiced one.
        if (bufferFilling > 0)
            bufferFilling = bufferLimit;
        end

        % Now we can resume standard unvoiced frame processing.
        if (inactiveFrameCounter < hangoverThreshold)
            % If we haven't padded the voiced frame with enough unvoiced
            % frames, do so by marking the frame voiced.
            iVad(frameCounter) = 1;
            inactiveFrameCounter = inactiveFrameCounter + 1;            
        elseif (inactiveFrameCounter == hangoverThreshold)
            % We've added the required unvoiced frame padding - set the
            % fade out flag for later processing
            applyFadeOut = 1;
            inactiveFrameCounter = inactiveFrameCounter + 1;
            firstVoicedFlag = 0;
        elseif (inactiveFrameCounter > hangoverThreshold)
            % Now that the padding and fade out are registered, mark the
            % end of unvoiced frame hangover smoothing so that we can add
            % silence in place of the noise.
            firstVoicedFlag = 0;
        end            
    else
        % Let's see if this is the first voice framed after a series of
        % unvoiced frames.
        if (firstVoicedFlag == 0)
            if (bufferFilling < bufferLimit)
                % If this is the first frame, we store it in the buffer,
                % because we only want to add voice sections that are a few
                % frames long (to compensate for anomalous spikes in
                % energy).
                bufferVoicedFrames = [bufferVoicedFrames; frame'];
                bufferFilling = bufferFilling + 1;
                
                if (bufferFilling == bufferLimit)
                    % We have filled the buffer, so this is a legit voiced
                    % frame! Add it to the output signal and clear the buffer.
                    % Also apply the fade in envelope and reset the first
                    % voiced flag.
                    frame = bufferVoicedFrames';
                    frame(1:frameLengthSamples) = frame(1:frameLengthSamples) .* fadeInEnvelope;
                    firstVoicedFlag = 1;
                    bufferFull = 1;
                end
            end
        else        
            % If this is simply the middle of an array of voiced frames, we
            % keep the unvoiced frame counter and fade out flags clean.
            inactiveFrameCounter = 0;
            applyFadeOut = 0;
        end
    end

    % Build threshold marker (debugging)
    thresholdMarker = [thresholdMarker; threshold * ones(length(frame), 1)];
    
    % Construct the VAD Marker and the final signal 
    if (iVad(frameCounter) == 1)
        if (bufferFilling > 0)
            if (bufferFull == 1)
                % The buffer is full of voiced frames that need to be used
                vadMarker = [vadMarker; 0.5 * ones(length(frame), 1)];
                outputSignal = [outputSignal; frame'];
                
                % Reset buffer
                bufferFilling = 0;
                bufferFull = 0;
                bufferVoicedFrames = [];
            end        
        else
            % Build VAD marker and output signal
            vadMarker = [vadMarker; 0.5 * ones(length(frame), 1)];
            outputSignal = [outputSignal; frame'];
        end
    else
        if (bufferFilling == bufferLimit)
            % Zero out the buffer contents and add them 
            frame = [bufferVoicedFrames; frame'];
            bufferFilling = 0;
            bufferVoicedFrames = [];
        end
        
        if (applyFadeOut == 1)
            % Fade out the previous unvoiced frames before adding silence
            fadeOutEnd = length(outputSignal);
            fadeOutStart = (fadeOutEnd - (frameLengthSamples * hangoverThreshold)) + 1;            
            outputSignal(fadeOutStart:fadeOutEnd) = outputSignal(fadeOutStart:fadeOutEnd) .* fadeOutEnvelope';
        end

        % Build VAD marker and output signal
        vadMarker = [vadMarker; zeros(length(frame), 1)];
        outputSignal = [outputSignal; zeros(length(frame), 1)];
    end
    
    % Scale Emin and Emax to compensate for anomalies in energy spectra
    energyMin = energyMin * deltaEmin;
    energyMax = energyMax * deltaEmax;
    
end

% Process output
outputSignal = processAudio(outputSignal);

% Write the final output voice-only wav to disk, and sound it
wavwrite(outputSignal, fs, 'VoiceOnlyOutput.wav');
% sound(outputSignal, fs);

% Plots

% Original signal (Voice signal + noise mix)
% Note: we exaggerate the noise in the plot a bit for increased visibility
if (nargin > 1)
    subplot(2, 1, 1), plot(noiseSignal * 2, 'k'), grid on, hold on, plot(voiceSignal, 'g'), plot(vadMarker, 'r', 'Linewidth', 2), grid on;
    title('Speech Signal with Voice Activity Detection'),xlabel('Time'), ylabel('Amplitude');
else
    subplot(2, 1, 1), plot(inputSignal, 'g'), hold on, plot(vadMarker, 'r', 'Linewidth', 2), grid on
    title('Speech Signal with Voice Activity Detection'),xlabel('Time'), ylabel('Amplitude');
end

% Reconstructed signal
hold off
subplot(2, 1, 2), plot(outputSignal), grid on, axis tight
title('Reconstructed Voice Only signal'),xlabel('Time'), ylabel('Amplitude');

end


% *************************************************************************
% Small helper function to process audio - mono conversion, DC removal, etc
% *************************************************************************
function in = processAudio(inp)

% Mono-fy, remove DC, normalize
inp = mean(inp, 2);
inp = inp - mean(inp);
inp = 0.99*inp/max(abs(inp));
in = inp;

end


% *************************************************************************
% Small helper function to resample audio to a higher rate
% *************************************************************************
function up = upSample(inp, origFs, targetFs)

ratio = targetFs/origFs;
ratioNumerator = round(ratio * 10000);
up = resample(inp, ratioNumerator, 10000);

end
