function [S0,a,b,theta] = postProc(S0,a,b,theta,TE,TR)
%[S0,a,b,theta] = postProc(S0,a,b,theta)
%Post processing step that can be used to get consistent estimates of the
%parameters. (Not needed for band removal)
%Input:
% S0, a, b, theta  -  Estimated parameter values
% TE, TR  -  Echo and repitition times.
%Output:
% S0, a, b, theta  -  Updated parameter estimates 
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

%False negative optimum
I1= b < 0 && b < 0;
  a(I1)=-a(I1);
  b(I1) = -b(I1);
  S0(I1)=S0(I1)*exp(1i*(-1)^(theta(I1)<0)*pi*TE/TR);
  theta(I1) = theta(I1) + (-1)^(theta(I1)>0)*pi; %We might be in either of the minima

% Wrapping problem
I2=theta>pi;
  k=round(theta(I2)/(2*pi));
  theta(I2)=theta(I2)-2*pi*k;
  S0(I2)=S0(I2)*exp(1i*2*pi*k*TE/TR);

I3=theta<-pi;
  k=round(theta(I3)/(-2*pi));
  theta(I3)=theta(I3)+2*pi*k;
  S0(I3)=S0(I3)*exp(-1i*2*pi*k*TE/TR);
