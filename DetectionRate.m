%{
Author : Mohammad Hovaidi Ardestani
Date and Place: 10.11.2015, CIN, Tuebingen, Germany
the last modification date : 29.01.2016
.
.
.
This function calculates a detection rate based on the percentage of time that 
visual activity is above a defined threshold
%}


function [Detections, DetectionRate] = DetectionRate(VsnOutput, WinSize, Threshold,SNR,delay)

% keyboard
DetectedSignal = VsnOutput > Threshold ;
% Signal = VsnOutput > 1e-2;
Detections = sum(sum(DetectedSignal));
DR = zeros(1,length(VsnOutput));
for i = 1:WinSize:length(VsnOutput)- WinSize
        if (all (any (DetectedSignal(:,i:i+WinSize) == 1)));
        DR(1,i:i+WinSize) = 1;
        end
end
plot(DR);
ylim([-0.1,1.1]);
DetectionRate = (sum(DR==1)/(length(DR)-delay))*100;

% DetectionRate = (Detections/sum(sum(Signal)))*100;
% title(['motion recognition rate[Baseline] = ',num2str(DetectionRate),'%',' with \bf(\sigma) = ',num2str(SNR)]);
title(['motion recognition rate[Synchronous] = ',num2str(DetectionRate),'%',' with \bf(\sigma) = ',num2str(SNR)]);
% title(['motion recognition rate[Asunchronous] = ',num2str(DetectionRate),'%',' with \bf(\sigma) = ',num2str(SNR)]);
end
