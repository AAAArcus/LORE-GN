function im_out=mat2image(im)
%im_out=mat2image(im)
%Converts a 3D array of images to one big image with all images put next to each
%other row-wise. The outputed image is adapted to the screen to maintain 
%scaling of the original image.
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

[y, x, z]=size(im);

screen=get(0,'screensize');

horizontal=round((0.7*screen(3))/x);
vertical=round((0.7*screen(4))/y);

nbPics=horizontal*vertical;

while nbPics<z
  horizontal=horizontal+1;
  vertical=vertical+1;
  nbPics=horizontal*vertical;
end
if z<=horizontal
  im_out=zeros(y,z*x);
else
  im_out=zeros(ceil(z/horizontal)*y,horizontal*x);
end

kk=1;
mm=1;
nn=1;
while kk<=z
  im_out(mm:(mm-1+y),nn:(nn-1+x))=im(:,:,kk);
  kk=kk+1;
  nn=nn+x;
  
  
  if nn==y*horizontal+1
    mm=mm+y;
    nn=1;
  end
end