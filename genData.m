function data=genData(S0,a,b,theta,deltaThetas,TE,TR)
%data=gendata(S0,a,b,theta,deltaThetas,TE,TR)
%Generate data for the bSSFP model and given parameters for each
%experiment phase increment in deltaThetas (radians). 
%Used for simulation and nonlinear least squares fitting.
%Input: 
% S0,a,b,theta  -  Model parameters
% deltaThetas  -  vector of image phase increments in radians
% TE, TR  -  Echo and repitition times
%Output:
% data  -  Generated data vector
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

data = S0.*exp(1i.*theta.*TE/TR).*(1-a.*exp(-1i.*(theta+deltaThetas)))./(1-b.*cos(theta+deltaThetas));

