function [thetaEst,S0Est,aEst,bEst,loss]=LMORM(data,deltaThetas,x0,TE,TR)
%[thetaEst,MEst,aEst,bEst,loss]=LMORM(data,deltaThetas,x0,TE,TR)
%Compute bSSFP model parameters using Levenberg-Marquardt (LM) based on 
%MATLABs built in lsqnonlin function.
%
%Input:
% data  -  Vector of complex-valued data to be fitted to.
% deltaThetas  -  Vector of corresponding phase increments (radians)
% x0 = [theta; real(S0); imag(S0); a; b]  -  Start guess for LM
% TE, TR  -  Echo and repitition times.
%Output:
% thetaEst,S0Est,aEst,bEst  -  Parameter estimates
% loss  -  criterion function value
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

% Model error and loss
error = @(x) [real(data-genData(x(2)+1i*x(3),x(4),x(5),x(1),deltaThetas,TE,TR));imag(data-genData(x(2)+1i*x(3),x(4),x(5),x(1),deltaThetas,TE,TR))];

%Levenberg-Marquardt implementation from MATLAB lsqnonlin function
OPTIONS = optimset('Algorithm','levenberg-marquardt','Display','off');
[x,loss]=lsqnonlin(error,x0,[],[],OPTIONS);

%Set output parameters
S0Est=x(2)+1i*x(3);
thetaEst=x(1);
bEst = x(5);
aEst=x(4);

end