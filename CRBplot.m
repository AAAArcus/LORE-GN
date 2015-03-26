%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the SNR needed to obtain a specified relative standard deviation
% (Pt) for T1 and T2. Generates Fig. 3 in the paper below.
%
% Please refer to the paper for more details: 
% "Parameter estimation approach to banding artifact reduction in balanced steady-state free precession"
% http://onlinelibrary.wiley.com/doi/10.1002/mrm.24986/full
%
% Written by: Marcus Björk, R. Reeve Ingle, Erik Gudmundson, Joëlle K. Barral, 2013 
%
% Copyright, Marcus Björk, Uppsala University, 2013
% 
% This file is part of the LORE-GN package.
% 
% LORE-GN is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% LORE-GN is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with LORE-GN.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

%User parameters
TR = 5; %ms
TE = 2.5;
alpha = 30*pi/180; %deg
KM0=exp(1i*0*pi/180);
theta = 0*pi/180;
deltaThetas=[0 pi/2 pi 3*pi/2]';
N1=1000; % %grid size in T1 direction. In paper: 2901
N2=100; %grid size in T2 direction. In paper: 196
Pt=0.05;  %Threshold for needed SNR, fraction of the true parameter value

%Compute constant parameters
ca = cos(alpha);
sa = sin(alpha);

%Create T1 and T2 grids
T1=linspace(100,3000,N1);
T2=linspace(5,200,N2);

h = waitbar(0, 'Computing, please wait...');

%Pre-allocation
signalMag=zeros(N1,N2);
CRBmat=zeros(N1,N2,4);
tic;
for k2=1:N2
    for k1=1:N1
        %Compute current parameters
        E1 = exp(-TR/T1(k1));
        E2 = exp(-TR/T2(k2));
        a = E2;
        b = E2*(1 - E1*ca - E1 + ca)/(1 - E1*ca - (E1-ca)*E2^2);
        KM = 1i*KM0*sa*( 1 - E1 )/( 1 - E1*ca - (E1-ca)*E2^2);
        %Generate data
        dataTmp=genData(KM*exp(-TE/T2(k2)),a,b,theta,deltaThetas,TE,TR);
        if T2(k2)>T1(k1) %Set nan for T2>T1 (not practical)
            CRBmat(k1,k2,:)=nan;
        else
            signalMag(k1,k2)=mean(abs(dataTmp));
            CRBmat(k1,k2,:)=sqrt(CRBT1T2(T1(k1),T2(k2),real(KM0),imag(KM0),theta,alpha,deltaThetas,TE,TR,1));
        end
    end
    waitbar(k2/N2, h, ['Computing, please wait... Time left: ' num2str(round((N2-k2)*toc/k2)) ' s']);
end
delete(h);

%MRI SNR
SNRneededT1Mat=signalMag.*CRBmat(:,:,1)./(Pt*repmat(T1',1,length(T2)));
SNRneededT2Mat=signalMag.*CRBmat(:,:,2)./(Pt*repmat(T2,length(T1),1));
SNRneededKM0Mat=signalMag.*CRBmat(:,:,3)./abs(Pt*KM0);
SNRneededThetaMat=signalMag.*CRBmat(:,:,4)./(Pt*pi/2);


%% MRI_SNR
%T1: needed SNR for relative std of Pt (see variable)
figure
imagesc(T2,T1,log10(SNRneededT1Mat),[1.4 2.6]);
set(gca,'YDir','normal','YTick',[100 500 1000 1500 2000 2500 3000],'XTick',[5 50 100 150 200]);
colormap(jet(256));
xlabel('T2 [ms]')
ylabel('T1 [ms]')
h=colorbar;
set(h,'ylim',[1.41 2.6],'YTick',[1.477 1.6 1.78 2 2.176 2.3 2.5],'YTickLabel',[30 40 60 100 150 200 300])
ylabel(h,'SNR (log-scale)')
set(gcf,'position',[500 500 400 296])

%T2: needed SNR for relative std of Pt (see variable)
figure
imagesc(T2,T1,log10(SNRneededT2Mat),[1.4 2.6]);
set(gca,'YDir','normal','YTick',[100 500 1000 1500 2000 2500 3000],'XTick',[5 50 100 150 200]);
colormap(jet(256));
xlabel('T2 [ms]')
ylabel('T1 [ms]')
h=colorbar;
set(h,'ylim',[1.41 2.6],'YTick',[1.477 1.6 1.78 2 2.176 2.3 2.5],'YTickLabel',[30 40 60 100 150 200 300])
ylabel(h,'SNR (log-scale)')
set(gcf,'position',[500 500 400 296])

