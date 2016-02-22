%{
Author : Mohammad Hovaidi Ardestani
Date and Place: 19.08.2015, CIN, Tuebingen, Germany
the last modification date : 29.01.2016
.
.
.
This function defines a filter through which the activity of visual field is 
being processed and passed to the motor field as the external stimulus
%}


function WW = VsnToMtrKernel(NN, gama, alpha, C)

% fac = 1e-3;  %%% default
fac = 100e-4;
%%% defines the strength of the filter.

pat = zeros(NN,NN);

    for location = 1:NN;  
     j = (1:NN)'; 
     dis = min(abs(j-location),NN-abs(j-location)); 
     pat(:,location)= ((1 + cos(alpha*dis))/2).^gama;
    end
    
WW = pat*pat'; WW = WW/WW(1,1); WW = fac *(WW-C);
% WW = 1e-3*(WW-C);
WW = circshift(WW, [0, 0]); %%% Direction From the last towards the first ensemble   
% plot(WW(:,NN/2));
% hold on
% P = WW(:,NN/2)';
% Q = 1:1:NN;
% idx = (-C <= P & P <= 0);
% hold on, plot(Q(idx), P(idx), '*r')
% hold off
% plot(pat)
end