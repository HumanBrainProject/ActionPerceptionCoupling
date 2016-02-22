%{
Author : Mohammad Hovaidi Ardestani
Date and Place: 03.07.2015, CIN, Tuebingen, Germany
the last modification date : 29.01.2016
.
.
.
This function implements a Gaussian filter to smooth the population activities 
of each ensemble.
%}


function GaussianFilter =  GaussianFilter(sigma,G_size)
    x = linspace(-G_size / 2, G_size / 2, G_size);
    GaussianFilter = exp(-x .^ 2 / (2 * sigma ^ 2));
    GaussianFilter = GaussianFilter / sum (GaussianFilter);
%   GaussianFilter = 1/sqrt(2*pi*sigma*2) * exp(-(x - mu).^2/(2*sigma^2));
%   plot(GaussianFilter)
end