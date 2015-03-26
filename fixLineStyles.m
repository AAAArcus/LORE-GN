function fixLineStyles(h)
%fixLineStyles(h)
%Fix the line-styles (set markers) for CRBvsSNR plots where the CRB is assumed to be
%plotted last (black solid). If no figure handle is given, gca is used. The 
%function has no output parameters.
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

if nargin==1
    lines=get(h,'children');
else
    lines=get(gca,'children');
end

N=length(lines);

markerCell={'x','+','o','.'};
markerSize=[6 6 8 6];
for k=2:N
    set(lines(k),'marker',markerCell{k-1},'markerSize',markerSize(k-1));%,'LineStyle',styleCell{k-1});
end
