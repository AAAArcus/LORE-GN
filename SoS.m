function imageSoS=SoS(im)
%imageSoS=SoS(im)
%Calculates the sum of squares image from a 3 dimensional matrix of images
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

[nbrow, nbcol, thetaLen]=size(im);

imageSoS=zeros(nbrow,nbcol);
for kk=1:thetaLen
  imageSoS=imageSoS+abs(im(:,:,kk)).^2;
end
imageSoS=1/sqrt(thetaLen)*sqrt(imageSoS);