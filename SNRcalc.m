function [SNR_lin, SNR_MRI, SNRim, noiseVar]=SNRcalc(im,mask)
%[SNR_lin, SNR_MRI, SNRim, noiseVar]=SNRcalc(im,mask)
%Calculate average SNR for an image matrix (possibly 3D) using both
%standard definition and MRI definition.
%Input:
% im  -  image matrix
% mask  -  mask for the background
%Output:
% SNR_lin  -  SNR in linear scale: root(mean(abs(signal)^2))/std
% SNR_MRI  -  SNR in the MRI definition: mean(abs(signal))/std
% SNRim  -  SNR image averaged over phase cycles only (matrix)
% noiseVar  -  Estimated background noise variance
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

[x,y,z]=size(im);

if nargin==1
    mask=ones(x,y);
end

%To avoid the transition region between tissue and background the masks are
%shrunken and grown by a few pixels pixel. Otherwise the SNR depend too
%much on the threshold. These parameters will be data dependent to some
%extent.
maskSignal=growShrink(mask,-5);
maskNoise=~growShrink(maskSignal,30);

signal=abs(im(repmat(maskSignal,[1 1 z])));
noise=im(repmat(maskNoise,[1 1 z]));
SNRim=zeros(x,y);
if isempty(noise)
    sigma=input('No noise standard deviation could be calculated, since there is no masked background. Input std: ');
    SNR_lin=sqrt(mean(signal.^2))/sigma;
    SNR_MRI=mean(signal)/sigma;
    noiseVar=sigma^2;
    SNRim(maskSignal)=mean(reshape(signal,[],z),2)./sigma;
else
    SNR_lin=sqrt(mean(signal.^2))/std(noise);
    SNR_MRI=mean(signal)/std(noise);
    noiseVar=var(noise);
    SNRim(maskSignal)=mean(reshape(signal,[],z),2)./std(noise);
end