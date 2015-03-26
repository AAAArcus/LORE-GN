function CRBvec=CRB(x,noiseVar,deltaThetas,TE,TR)
%CRBvec=CRB(x,noiseVar,deltaThetas,TE,TR)
%Computes the numerical Cramér-Rao lower bound in terms of variance of the parameters.
%Input:
% x = [theta real(S0) imag(S0) a b]  -  the vector of true parameter values
% noiseVar  -  the noise variance 
% deltaThetas  -  vector of image phase increments in radians 
% TE, TR  -  Echo and repitition times
%Output:
% CRBvec = [theta_var S0r_var S0i_var a_var b_var]  -  Vector of minimum variances 
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

N=length(deltaThetas);
%Pre-allocation
FIM=zeros(5,5);

for k=1:N
  d=deltaThetas(k);
  %Derivatives
  dS_dtheta=- (x(4)*exp((TE*x(1)*1i)/TR)*(x(2) + x(3)*1i)*1i)/(exp(d*1i + x(1)*1i)*(x(5)*cos(d + x(1)) - 1)) + (TE*exp((TE*x(1)*1i)/TR)*(x(4)/exp(d*1i + x(1)*1i) - 1)*(x(2) + x(3)*1i)*1i)/(TR*(x(5)*cos(d + x(1)) - 1)) + (x(5)*exp((TE*x(1)*1i)/TR)*sin(d + x(1))*(x(4)/exp(d*1i + x(1)*1i) - 1)*(x(2) + x(3)*1i))/(x(5)*cos(d + x(1)) - 1)^2;
  dS_dS0r=(exp((TE*x(1)*1i)/TR)*(x(4)/exp(d*1i + x(1)*1i) - 1))/(x(5)*cos(d + x(1)) - 1);
  dS_dS0i=(exp((TE*x(1)*1i)/TR)*(x(4)/exp(d*1i + x(1)*1i) - 1)*1i)/(x(5)*cos(d + x(1)) - 1);
  dS_da=(exp((TE*x(1)*1i)/TR)*(x(2) + x(3)*1i))/(exp(d*1i + x(1)*1i)*(x(5)*cos(d + x(1)) - 1));
  dS_db=-(exp((TE*x(1)*1i)/TR)*cos(d + x(1))*(x(4)/exp(d*1i + x(1)*1i) - 1)*(x(2) + x(3)*1i))/(x(5)*cos(d + x(1)) - 1)^2;
  
  %Jacobian
  J=[dS_dtheta;dS_dS0r;dS_dS0i;dS_da;dS_db];
  
  %Update normalized Fisher Information Matrix (FIM)
  FIM=FIM+real(J*J');
end
%Scale by noise variance
FIM=2/noiseVar*FIM;

%Return diagonal only (no covariance)
CRBvec=diag(inv(FIM));