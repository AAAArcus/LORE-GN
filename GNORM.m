function [thetaEst,S0Est,aEst,bEst,loss,ii] = GNORM(data,deltaThetas,x0,TE,TR,tol,maxIter)
%[thetaEst,S0Est,aEst,bEst,loss,ii] = GNORM(data,deltaThetas,x0,TE,TR,tol,maxIter)
%Gauss Newton (GN) algorithm for minimizing the NLS criterion of the bSSFP model.
%
%Input:
% data  -  Vector of complex-valued data to be fitted to.
% deltaThetas  -  Vector of corresponding phase increments (radians)
% x0 = [theta; real(S0); imag(S0); a; b]  -  Start guess for GN
% tol  -  tolerance  of the norm of the gradient (maximum allowable)
% maxIter  -  Maximum number of iterations.
% TE, TR  -  Echo and repitition times.
%Output:
% thetaEst, S0Est, aEst, bEst  -  Parameter estimates
% loss  -  Vector of criterion function values over iterations
% ii  -  The number of iterations made
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

%Initialize
%warning('off','MATLAB:rankDeficientMatrix'); %Turns off warnings for rank deficiency if GN is applied to low SNR data
y = [real(data); imag(data)];
nbrData=length(data);
x_new = x0(:);
loss=zeros(maxIter+1,1);
normGrad=tol;

ii=1;
while normGrad>=tol && ii<=maxIter
  x=x_new;
  %Compute the Jacobian and crit function with current estimates x
  theta_n = x(1) + deltaThetas;

  %Model and criterion (loss) value at current step
  f_kTmp = genData((x(2)+1i*x(3)),x(4),x(5),x(1),deltaThetas,TE,TR);
  f_k = [real(f_kTmp); imag(f_kTmp)];
  loss(ii)=1/nbrData*sum(abs(y-f_k).^2);    %Current loss

  %Common factors in the jacobian (for speed and simplified notation)
  c_n = cos(theta_n);
  s_n = sin(theta_n);
  c2_n=cos(theta_n-x(1)*TE/TR);
  s2_n=sin(theta_n-x(1)*TE/TR);
  
  %Jacobian elements
  dfre_n_dx1 = (c2_n.*(x(3)*c_n + x(2)*s_n) - s2_n.*(x(2)*c_n - x(3)*s_n) + (c2_n*(TE - TR).*(x(3)*c_n - x(3)*x(4) + x(2)*s_n))/TR + (s2_n*(TE - TR).*(x(2)*x(4) - x(2)*c_n + x(3)*s_n))/TR)./(x(5)*c_n - 1) + (x(5)*s_n.*(c2_n.*(x(2)*x(4) - x(2)*c_n + x(3)*s_n) - s2_n.*(x(3)*c_n - x(3)*x(4) + x(2)*s_n)))./(x(5)*c_n - 1).^2;
  dfre_n_dx2 = -(s2_n.*s_n - c2_n.*(x(4) - c_n))./(x(5)*c_n - 1);
  dfre_n_dx3 = (sin((TE*x(1))/TR) + x(4)*s2_n)./(x(5)*c_n - 1);
  dfre_n_dx4 = (x(2)*c2_n + x(3)*s2_n)./(x(5)*c_n - 1);
  dfre_n_dx5 = -(c_n.*(c2_n.*(x(2)*x(4) - x(2)*c_n + x(3)*s_n) - s2_n.*(x(3)*c_n - x(3)*x(4) + x(2)*s_n)))./(x(5)*c_n - 1).^2;
  dfim_n_dx1 = - (c2_n.*(x(2)*c_n - x(3)*s_n) + s2_n.*(x(3)*c_n + x(2)*s_n) - (c2_n*(TE - TR).*(x(2)*x(4) - x(2)*c_n + x(3)*s_n))/TR + (s2_n*(TE - TR).*(x(3)*c_n - x(3)*x(4) + x(2)*s_n))/TR)./(x(5)*c_n - 1) - (x(5)*s_n.*(c2_n.*(x(3)*c_n - x(3)*x(4) + x(2)*s_n) + s2_n.*(x(2)*x(4) - x(2)*c_n + x(3)*s_n)))./(x(5)*c_n - 1).^2;
  dfim_n_dx2 = -(sin((TE*x(1))/TR) + x(4)*s2_n)./(x(5)*c_n - 1);
  dfim_n_dx3 = -(s2_n.*s_n - c2_n.*(x(4) - c_n))./(x(5)*c_n - 1);
  dfim_n_dx4 = (x(3)*c2_n - x(2)*s2_n)./(x(5)*c_n - 1);
  dfim_n_dx5 = (c_n.*(c2_n.*(x(3)*c_n - x(3)*x(4) + x(2)*s_n) + s2_n.*(x(2)*x(4) - x(2)*c_n + x(3)*s_n)))./(x(5)*c_n - 1).^2;
  
  %The upper half is constructed from the derivatives of the real part,
  %the lower half from the imaginary part of the data model
  Jf = [dfre_n_dx1 dfre_n_dx2 dfre_n_dx3 dfre_n_dx4 dfre_n_dx5   ;...
    dfim_n_dx1 dfim_n_dx2 dfim_n_dx3 dfim_n_dx4 dfim_n_dx5];

  %Take one step in the G-N equation
  %Backtracking steplength
  alpha=2;
  okStep=0;
  mu=0.1;     
  while ~okStep
    alpha=alpha/2;
    p=(Jf\(y - f_k));
    x_new = x +  alpha*p;
    f_knewTmp=genData((x_new(2)+1i*x_new(3)),x_new(4),x_new(5),x_new(1),deltaThetas,TE,TR);
    f_knew=[real(f_knewTmp);imag(f_knewTmp)];
    okStep=(sum((y-f_knew).^2)<=sum((y-f_k).^2)-mu*alpha*p'*Jf'*(y-f_k));  %f to be minimized?!
  end

  normGrad=norm(Jf'*(y-f_knew))/norm(y);   %norm of gradient for stopping condition
  ii=ii+1;
  
end
loss(ii)=1/nbrData*sum(abs(y-f_knew).^2);   %Final loss
loss=loss(1:ii);


%% Output

if loss(1)<loss(ii) %If G-N fails to decrease the criterion (instability), use initial estimate (LORE)
  x=x0;
end
thetaEst = x(1);
S0Est = x(2)+1i*x(3);
aEst = x(4);
bEst = x(5);


