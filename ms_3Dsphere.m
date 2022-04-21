function [ h ] = ms_3Dsphere( r, sc)
%MS_3DGAUSSIAN Summary of this function goes here
%   Detailed explanation goes here
if nargin<1; r=3; end
if nargin<2; sc=[1 1 1]; end
N = 2*r+1; %// Define size of Gaussian mask

%// Generate Gaussian mask
ind = -floor(N/2) : floor(N/2);
[X, Y, Z] = meshgrid(ind, ind, ind); 
X = X*sc(1); Y = Y*sc(2); Z = Z*sc(3);
% h = (2*pi)^(-3/2) * sigma^-3 * exp(-(X.^2 + Y.^2 + Z.^2) / (2*sigma*sigma));
h = (X.^2 + Y.^2 + Z.^2);

h(h<=r^2)=1; h(h>r^2)=0; 

end

