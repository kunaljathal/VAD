# VAD
Voice Activity Detection System

VoiceActivityDetector implements a Vocal Activity Detection algorithm. 
It works by using voiced band energy ratios, periodicity measures, running 
and dynamic minimum and maximum RMS energy estimates, adaptive and noise 
resistant threshold computation, and hangover smoothing to fade in and out 
speech/noise boundaries. The final output is a reconstructed voice only audio 
file, together with a visual representation showing the VAD marked sections 
and the voice plot

Usage

VoiceActivityDetector(input, noise) - Performs VAD on the input audio
signal 'input'. 'noise' is an optional wav of noise to mix into the
input audio file.

For input, the following voice audio files are present in the folder:

male1.wav     -       sample male voice
male2.wav     -       sample male voice
male3.wav     -       sample male voice
female1.wav   -       sample female voice
female2.wav   -       sample female voice

For noise, the following files are present in the folder:

'white.wav'   -       standard white noise, constant level
'pink.wav'    -       standard pink noise, constant level
'babble.wav'  -       generic people talking, dynamically changing levels
'music.wav'   -       music used as noise, dynamically changing levels

The noise file is optional. You can simply provide your own already noisy
signal directly as 'input' if you'd like, or you can also provide your own
individual voice and noise files separately.

NOTE: If you are using only a single noisy input file, please make
sure the noise is reasonably softer than the voice input. The VAD works
best when the noise is around 10of the signal amplitude. Also, ensure
the first ~300ms of the audio is only noise (no voice).

Finally , to play around with the accuracy of the VAD, you can try tweaking
the following variables directly in the code itself:

frameLength           -   ~10 to 20ms (if you want to go below 10, you 
                                       will need to change the upperLagMS 
                                       value in Periodicity.m to be lower)
                                         
hangOverThreshold     -   ~3 to 10    (# of unvoiced frames to pad  
                                       unconditionally after the end of a 
                                       voiced frame

noiseScalingFactor    -   ~0.05 to ~0.2 (amplification of noise compared
                                         to voice signal expressed as a 
                                         percentage

