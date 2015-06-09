% *************************************************************************
% Kunal Jathal
% N19194426
% DST 2 - Final Project
% 
% Name:     frameSplit
%
% Purpose:  Helper Function
%
% Description:
%
% This splits the vector inputSignal(:) up into frames. Each frame is of 
% length LEN and occupies one row of the output matrix. The last few frames
% of X will be ignored if its length is not divisible by LEN. It is an error 
% if X is shorter than LEN.
%
%
% Usage
% 
% F = frameSplit(input, windowLength, increment) has frames beginning at 
% increments of 'increment'. The centre of frame I is 
% input((I-1)*increment+(windowLength+1)/2) for I=1,2,...
% The number of frames is fix((length(input)-windowLength+increment)/increment)
%
% F = frameSplit(input, windowLength) or 
% F = frameSplit(input, windowLength, increment) multiplies
% each frame by WINDOW(:)
% 
% Note: This code was inspired by the enframe function, which is open
% source and part of the voicebox toolbox
% 
% *************************************************************************

function finalSignal = frameSplit(input, windowLength, increment)
 
% Read input variables
inputSize=length(input(:));
winLength=length(windowLength);

% Initialize window length
if (winLength == 1)
    len = windowLength;
else
    len = winLength;
end

% Number of arguments check
if (nargin < 3)
   increment = len;
end

% Window processing! Figure out the frame segmentation etc.
numFrames = fix((inputSize-len+increment)/increment);
finalSignal=zeros(numFrames,len);
indexFactor= increment*(0:(numFrames-1)).';
indices = (1:len);
finalSignal(:) = input(indexFactor(:,ones(1,len))+indices(ones(numFrames,1),:));

% Construct final signal
if (winLength > 1)
    winFinal = windowLength(:)';
    finalSignal = finalSignal .* winFinal(ones(numFrames,1),:);
end

return