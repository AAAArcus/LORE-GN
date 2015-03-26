function [thetaEst,S0Est,aEst,bEst] = LORE(data,deltaThetas,TE,TR)
%[thetaEst,S0Est,aEst,bEst] = LORE(data,deltaThetas,TE,TR)
%Find bSSFP model parameters by the use of LS for all parameters by using a data 
%dependent regressor matrix and a overparameterized model (6 parameters).
%
%Input:
% data  -  Vector of complex-valued data to be fitted to.
% deltaThetas  -  Vector of corresponding phase increments (radians)
% x0 = [theta; real(S0); imag(S0); a; b]  -  Start guess for GN
% TE, TR  -  Echo and repitition times.
%Output:
% thetaEst,S0Est,aEst,bEst  -  Parameter estimates
%
%Written by: Marcus Björk, R. Reeve Ingle, Erik Gudmundson, Joëlle K. Barral, 2013 
%
%Copyright, Marcus Björk, Uppsala University, 2013
% 
%This file is part of the LORE-GN package.
% 
%LORE-GN is free software: you can redistribute it and/or modify
%it under the terms of the GNU General Public License as published by
%the Free Software Foundation, either version 3 of the License, or
%(at your option) any later version.
% 
%LORE-GN is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
% 
%You should have received a copy of the GNU General Public License
%along with LORE-GN.  If not, see <http://www.gnu.org/licenses/>.

%Extract real and imaginary part of data
Ire = real(data);
Iim = imag(data);
dataLen = length(deltaThetas);

if (length(data) ~= dataLen)
  error('ERROR: Size of data and deltaThetas must be equal!')
end

%Model: yn=(alpha-beta*e^(i*deltaOmega*TR))/(1-Re{gamma*e^(i*deltaOmega*TR))
% Calculate regressor matrix
A = [ones(dataLen,1) zeros(dataLen,1) -cos(deltaThetas) -sin(deltaThetas) Ire.*cos(deltaThetas) -Ire.*sin(deltaThetas)
     zeros(dataLen,1) ones(dataLen,1) sin(deltaThetas) -cos(deltaThetas) Iim.*cos(deltaThetas) -Iim.*sin(deltaThetas) ];

%Parameter vector: x=[gamma_r gamma_i alpha_r alpha_i beta_r beta_i]';
%Least squares
x = A\[Ire;Iim];

%Compute 
alpha=x(1)+1i*x(2);
beta=x(3)+1i*x(4);
gamma=x(5)+1i*x(6);

%Non unique way of obtaining estimates of the original parameters
thetaEst=-angle(beta/alpha);
S0Est=alpha*exp(-1i*thetaEst*TE/TR);
aEst=abs(beta/alpha);
bEst=abs(gamma);   

