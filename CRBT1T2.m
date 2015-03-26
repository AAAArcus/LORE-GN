function CRBvec=CRBT1T2(T1,T2,KM0r,KM0i,theta,alpha,deltaThetas,TE,TR,noiseVar)
% CRBvec=CRBT1T2(T1,T2,KM0r,KM0i,theta,alpha,deltaThetas,TE,TR,noiseVar)
%Computes the numerical Cramér-Rao lower bound in terms of variance of the parameters
%for the T1 T2 based model.
%Input:
% T1,T2,KM0r,KM0i,theta  -  true parameter values
% alpha  -  flip angle in radians
% noiseVar  -  the noise variance 
% deltaThetas  -  vector of image phase increments in radians, 
% TE, TR  -  Echo and repitition times.
%Output:
% CRBvec = [T1_var T2_var KM0_var theta_var]  -  Vector of minimum variances 
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

%Constant parameters
ca=cos(alpha);
sa=sin(alpha);
deltaThetas=deltaThetas(:);

%Compute gradient    
dS_dT1= -(TR.*exp(- deltaThetas.*1i - theta.*1i).*exp((TE.*theta.*1i)./TR).*exp(-TE./T2).*exp(TR./T1).*exp(TR./T2).*sa.*(ca - 1).*(exp((2.*TR)./T2) - 1).*(KM0i - KM0r.*1i).*(exp(deltaThetas.*1i + theta.*1i).*exp(TR./T2) - 1))./(T1.^2.*(exp(TR./T2).*cos(deltaThetas + theta) + exp(TR./T1).*ca - exp((2.*TR)./T2).*ca + exp(TR./T1).*exp((2.*TR)./T2) + exp(TR./T2).*cos(deltaThetas + theta).*ca - exp(TR./T1).*exp(TR./T2).*cos(deltaThetas + theta) - exp(TR./T1).*exp(TR./T2).*cos(deltaThetas + theta).*ca - 1).^2);
dS_dT2= (exp(- deltaThetas.*1i - theta.*1i).*exp((TE.*theta.*1i)./TR).*exp(-TE./T2).*exp(TR./T2).*sa.*(exp(TR./T1) - 1).*(KM0r + KM0i.*1i).*(TE - TR - TE.*exp(TR./T1).*ca + TE.*exp((2.*TR)./T2).*ca + TR.*exp(TR./T1).*ca + TR.*exp((2.*TR)./T2).*ca - TE.*exp(TR./T1).*exp((2.*TR)./T2) - TR.*exp(TR./T1).*exp((2.*TR)./T2) - TE.*exp(deltaThetas.*1i + theta.*1i).*exp(TR./T2) + 2.*TR.*exp(deltaThetas.*1i + theta.*1i).*exp(TR./T2) - TE.*exp(deltaThetas.*1i + theta.*1i).*exp((3.*TR)./T2).*ca + TE.*exp(deltaThetas.*1i + theta.*1i).*exp(TR./T1).*exp((3.*TR)./T2) + TE.*exp(deltaThetas.*1i + theta.*1i).*exp(TR./T1).*exp(TR./T2).*ca - 2.*TR.*exp(deltaThetas.*1i + theta.*1i).*exp(TR./T1).*exp(TR./T2).*ca).*1i - exp(- deltaThetas.*1i - theta.*1i).*exp((TE.*theta.*1i)./TR).*exp(-TE./T2).*exp(TR./T2).*cos(deltaThetas + theta).*sa.*(exp(TR./T1) - 1).*(KM0r + KM0i.*1i).*(TE.*exp(TR./T2) + TE.*exp(TR./T2).*ca - TE.*exp(TR./T1).*exp(TR./T2) - TE.*exp(deltaThetas.*1i + theta.*1i).*exp((2.*TR)./T2) + TR.*exp(deltaThetas.*1i + theta.*1i).*exp((2.*TR)./T2) - TE.*exp(TR./T1).*exp(TR./T2).*ca - TE.*exp(deltaThetas.*1i + theta.*1i).*exp((2.*TR)./T2).*ca + TR.*exp(deltaThetas.*1i + theta.*1i).*exp((2.*TR)./T2).*ca + TE.*exp(deltaThetas.*1i + theta.*1i).*exp(TR./T1).*exp((2.*TR)./T2) - TR.*exp(deltaThetas.*1i + theta.*1i).*exp(TR./T1).*exp((2.*TR)./T2) + TE.*exp(deltaThetas.*1i + theta.*1i).*exp(TR./T1).*exp((2.*TR)./T2).*ca - TR.*exp(deltaThetas.*1i + theta.*1i).*exp(TR./T1).*exp((2.*TR)./T2).*ca).*1i)./(T2.^2.*(exp(TR./T2).*cos(deltaThetas + theta) + exp(TR./T1).*ca - exp((2.*TR)./T2).*ca + exp(TR./T1).*exp((2.*TR)./T2) + exp(TR./T2).*cos(deltaThetas + theta).*ca - exp(TR./T1).*exp(TR./T2).*cos(deltaThetas + theta) - exp(TR./T1).*exp(TR./T2).*cos(deltaThetas + theta).*ca - 1).^2);
dS_dKM0r= -(exp((TE.*theta.*1i)./TR).*exp(-TE./T2).*sa.*(exp(-TR./T1) - 1).*(exp(- deltaThetas.*1i - theta.*1i).*exp(-TR./T2) - 1).*1i)./(((exp(-TR./T2).*cos(deltaThetas + theta).*(ca - exp(-TR./T1) - exp(-TR./T1).*ca + 1))./(exp(-(2.*TR)./T2).*(ca - exp(-TR./T1)) - exp(-TR./T1).*ca + 1) - 1).*(exp(-(2.*TR)./T2).*(ca - exp(-TR./T1)) - exp(-TR./T1).*ca + 1));
dS_dKM0i= (exp((TE.*theta.*1i)./TR).*exp(-TE./T2).*sa.*(exp(-TR./T1) - 1).*(exp(- deltaThetas.*1i - theta.*1i).*exp(-TR./T2) - 1))./(((exp(-TR./T2).*cos(deltaThetas + theta).*(ca - exp(-TR./T1) - exp(-TR./T1).*ca + 1))./(exp(-(2.*TR)./T2).*(ca - exp(-TR./T1)) - exp(-TR./T1).*ca + 1) - 1).*(exp(-(2.*TR)./T2).*(ca - exp(-TR./T1)) - exp(-TR./T1).*ca + 1));
dS_dtheta= - (exp(- deltaThetas.*1i - theta.*1i).*exp((TE.*theta.*1i)./TR).*exp(-TE./T2).*exp(-TR./T2).*sa.*(exp(-TR./T1) - 1).*(KM0r + KM0i.*1i))./(((exp(-TR./T2).*cos(deltaThetas + theta).*(ca - exp(-TR./T1) - exp(-TR./T1).*ca + 1))./(exp(-(2.*TR)./T2).*(ca - exp(-TR./T1)) - exp(-TR./T1).*ca + 1) - 1).*(exp(-(2.*TR)./T2).*(ca - exp(-TR./T1)) - exp(-TR./T1).*ca + 1)) + (TE.*exp((TE.*theta.*1i)./TR).*exp(-TE./T2).*sa.*(exp(-TR./T1) - 1).*(KM0r + KM0i.*1i).*(exp(- deltaThetas.*1i - theta.*1i).*exp(-TR./T2) - 1))./(TR.*((exp(-TR./T2).*cos(deltaThetas + theta).*(ca - exp(-TR./T1) - exp(-TR./T1).*ca + 1))./(exp(-(2.*TR)./T2).*(ca - exp(-TR./T1)) - exp(-TR./T1).*ca + 1) - 1).*(exp(-(2.*TR)./T2).*(ca - exp(-TR./T1)) - exp(-TR./T1).*ca + 1)) - (exp((TE.*theta.*1i)./TR).*exp(-TE./T2).*exp(-TR./T2).*sin(deltaThetas + theta).*sa.*(exp(-TR./T1) - 1).*(KM0r + KM0i.*1i).*(exp(- deltaThetas.*1i - theta.*1i).*exp(-TR./T2) - 1).*(ca - exp(-TR./T1) - exp(-TR./T1).*ca + 1).*1i)./(((exp(-TR./T2).*cos(deltaThetas + theta).*(ca - exp(-TR./T1) - exp(-TR./T1).*ca + 1))./(exp(-(2.*TR)./T2).*(ca - exp(-TR./T1)) - exp(-TR./T1).*ca + 1) - 1).^2.*(exp(-(2.*TR)./T2).*(ca - exp(-TR./T1)) - exp(-TR./T1).*ca + 1).^2);

%Jacobian
J=[dS_dT1 dS_dT2 dS_dKM0r dS_dKM0i dS_dtheta];

%Compute CRB matrix
CRB=1/2*inv(real(J'*J));
%Return diagonal elements (no covariance)
CRBvec=[CRB(1,1) CRB(2,2) CRB(3,3)+CRB(4,4) CRB(5,5)]'*noiseVar(:)';

