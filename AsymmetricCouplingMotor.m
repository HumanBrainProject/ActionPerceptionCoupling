%{
Author : Mohammad Hovaidi Ardestani
Date and Place: 07.08.2015, CIN, Tuebingen, Germany
the last modification date : 29.01.2016
.
.
.
This Function takes the number of Ensembles (NN),Gama,alpha, and C and then 
returns the Coupling Matrix of the field.
%%% The Parameter "C" determines the Iceberg Peak! %%%
%}


function WW = AsymmetricCouplingMotor(NN, gama, alpha, C)
fac = 200e-9; 
% in case of traveling stimulus this factor should be decreased to 20e-9.
pat = zeros(NN,NN);

    for location = 1:NN;  
     j = (1:NN)'; 
     dis = min(abs(j-location),NN-abs(j-location)); 
     pat(:,location)= ((1 + cos(alpha*dis))/2).^gama;
    end
    
WW = pat*pat'; 
WW = WW/WW(1,1); 
WW = fac*(WW-C); %%default
WW = circshift(WW, [1, 0]); %%% Direction From the last towards the first ensemble   
% plot(WW(:,NN/2));
% hold on
% P = WW(:,NN/2)';
% Q = 1:1:NN;
% idx = (-C <= P & P <= 0);
% hold on, plot(Q(idx), P(idx), '*r')
% hold off
% plot(pat)
end