%{
Author : Mohammad Hovaidi Ardestani
Date and Place: 10.08.2015, CIN, Tuebingen, Germany
the last modification date : 29.01.2016
.
.
.
This function returns a travelling pulse used as the external motor stimulus.
(GO SIGNAL) 
%}


function GaussianInput = MtrInput(M_size,sigma1,ST_v)
dt = 1e-3;                             % time step
Time = 2;                              % simulation time
T = 0:dt:Time;                         % time step vector
% ST_m = ST_v + delay;                   % Defines the time when stimulus shown
ST_m = ST_v;
v = 6;                                 % Defines the wave speed

GaussianInput = zeros(M_size,length(T));

        xx1 = linspace(-M_size/2, M_size/2 ,M_size);
        GaussianInput(:,ST_m) =  exp(-xx1 .^ 2 / (2 * sigma1 ^ 2)); 
        GaussianInput(:,ST_m) = GaussianInput(:,ST_m) / sum (GaussianInput(:,ST_m));    

    for t = ST_m+1:v:length(T)

        GaussianInput(:,t) = circshift(GaussianInput(:,t-1),[1,0]);%%% OppositeDirection: From the first towards the last ensemble
%         GuassianInput(:,t) = circshift(GuassianInput(:,t-1),[29,0]);%%% SameDirection: From the last towards the first ensemble   
        for i = 0:v
            GaussianInput(:,t+i+1) = GaussianInput(:,t+i);
        end

    end
    GaussianInput(:,ST_m+10:end) = 0;  % This makes the stimulus transient. It sholud be commented out if a travelling pulse is needed   
end
