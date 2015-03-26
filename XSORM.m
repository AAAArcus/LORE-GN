function S0Est = XSORM(data)
%S0Est = XSORM(data)
%This method used the cross-solution method presented at ISMRM 2010 by
%by Xiang & Hoff to estimate S0 from data given by:
%     I_n = M*(1 - a*exp(i*theta_n)/(1 - b*cos(theta_n),
%where theta_n = theta + deltaTheta_n. theta is the off-resonance
%and deltaThetas = [0, pi/2, pi, 3pi/2]. The method estimates S0 as the 
%center of an ellipsoid.
%Input:
% data  -  4x1 vector with data (also works for 3D matrix)
% S0Est  -  estimated S0
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

I = data;
x = real(data);
y = imag(data);
%Make sure data consists of four image points
[~, ~, c]=size(data);
if c==1 
    if (length(data) ~= 4)
        error('ERROR: This method requires exactly 4 data points!')
    else %we have a vector och 4 datapoints
        S0Est = ( (x(1).*y(3) - x(3).*y(1)).*(I(2) - I(4)) ...
            - (x(2).*y(4) - x(4).*y(2)).*(I(1) - I(3))  ) ./ ...
            (x(1).*y(2) + x(2).*y(3) + x(3).*y(4) + x(4).*y(1) ...
            - x(1).*y(4) - x(4).*y(3) - x(3).*y(2) - x(2).*y(1) );
    end
elseif c~=4
    error('ERROR: This method requires exactly 4 images!')
else %we have a 3D matrix of 4 images
     S0Est = ( (x(:,:,1).*y(:,:,3) - x(:,:,3).*y(:,:,1)).*(I(:,:,2) - I(:,:,4)) ...
            - (x(:,:,2).*y(:,:,4) - x(:,:,4).*y(:,:,2)).*(I(:,:,1) - I(:,:,3))  ) ./ ...
            (x(:,:,1).*y(:,:,2) + x(:,:,2).*y(:,:,3) + x(:,:,3).*y(:,:,4) + x(:,:,4).*y(:,:,1) ...
            - x(:,:,1).*y(:,:,4) - x(:,:,4).*y(:,:,3) - x(:,:,3).*y(:,:,2) - x(:,:,2).*y(:,:,1) );
end