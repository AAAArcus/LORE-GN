function x0=startGuess(TE,TR,T1,T2,alpha,S0)
%x0=startGuess(TE,TR,T1,T2,alpha,S0)
%Generate a start guess of the bSSFP model parameters.
%Input
% T1, T2, S0  -  Guess based on tissue (expected values) 
% TE, TR  -  Echo and repitition times.
% alpha  -  Flip angle in radians.
%Output:
% x0 = [theta;real(S0);imag(S0);a;b]  -  Parameters for the reparameterized model
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

E1 = exp(-TR/T1);
E2 = exp(-TR/T2);
ca = cos(alpha);
a = E2;
b = E2*(1 - E1*ca - E1 + ca)/(1 - E1*ca - (E1-ca)*E2^2);
%S0 = 1i*M0*sin(alpha)*( 1 - E1 )/( 1 - E1*ca - (E1-ca)*E2^2);  %Alternative if M0 is given as input
theta=3*pi; %Any theta>pi will cause problems for LM
x0=[theta real(S0) imag(S0) a b].';