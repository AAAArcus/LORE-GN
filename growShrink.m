function mask=growShrink(mask,n)
%Grow or shrink the mask. The mask is enlarged by n pixels. If n is
%negative the mask in shrunk. A 1 indicate the desired and a 0 the masked.
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

[Nrow, Ncol]=size(mask);

if nargin==1 %Default value is grow by 1.
    n=1;
end
if n==0; %Nothing to be done
    return
else
    S=sign(n);
    P=n>0;
    for k=1:abs(n) %Recursive use of mask for growing/shrinking several times
        %% Across columns
        %colDiff indicates mask borders 
        colDiff=diff(mask);
        %Partial new mask
        mask([colDiff==S;false(1,Ncol)]|[false(1,Ncol);colDiff==-S])=P;
        
        %% Across rows
        %rowDiff indicates mask borders 
        rowDiff=diff(mask,[],2);
        %Full new mask
        mask([rowDiff==S false(Nrow,1)]|[false(Nrow,1) rowDiff==-S])=P;
    end
end