function imMI=maxIntensity(data,complex)
%imMI=maxIntensity(data,complex)
%Takes any 3D matrix of images and returns the maximum intensity image.
%Set complex=1 to return the complex image.
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

[imMI,i]=max(abs(data),[],3);

if nargin==2 && complex==1
    [nRow,nCol]=size(data);
    for k1=1:nRow
        for k2=1:nCol
            imMI(k1,k2)=data(k1,k2,i(k1,k2));
        end
    end
end
    