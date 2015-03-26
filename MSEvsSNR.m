%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run Monte Carlo simulations and compute rMSE. Compare to CRB to determine 
% the efficiency of the algorithms at different SNRs. 
% It is possible to change the true parameters simulated,
% the number of simulations, and what methods to include.
% The computations can be quite timeconsuming for a large number of
% simulations and methods.
%
%Please refer to the paper for more details: 
%"Parameter estimation approach to banding artifact reduction in balanced steady-state free precession"
%http://onlinelibrary.wiley.com/doi/10.1002/mrm.24986/full
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
%close all
%% Options
MC = 200;  %number of simulations (10000 used in the paper)
N_SNR = 10; %Number of SNRs sampled (30 used in the paper)
methodsToTest={
  'LORE'
  'LORE-GN'
  'LM'
  'LMpost'
};

plotResults=1;

%% True parameter values
TR = 5; %ms
TE = 2.5;
T1 = 675;
T2 = 75;
alpha = 30*pi/180; %deg
KM0=1;
theta = 0*pi/180;

T1start = 750;
T2start = 60;

E1 = exp(-TR/T1);
E2 = exp(-TR/T2);
ca = cos(alpha);
sa = sin(alpha);
a = E2;
b = E2*(1 - E1*ca - E1 + ca)/(1 - E1*ca - (E1-ca)*E2^2);
S0 = 1i*KM0*sa*( 1 - E1 )/( 1 - E1*ca - (E1-ca)*E2^2)*exp(-TE/T2);

%% Generate noise-free data, compute SNR
step=pi/2;
deltaThetas=(0:step:(2*pi-step))';
theta_n = theta + deltaThetas; %in radians
dataNoNoise=genData(S0,a,b,theta,deltaThetas,TE,TR);

%MRI SNR
signalMag=mean(abs(dataNoNoise));  %Power
SNR_MRI=10.^linspace(log10(5),log10(100),N_SNR);
noiseVar=(signalMag./SNR_MRI).^2;

%% Pre-allocation and initialization
startTimeTot = cputime;
noiseLen = length(noiseVar);
nbMethods=length(methodsToTest);
thetaLen = length(deltaThetas);
S0Est = zeros(noiseLen,MC,nbMethods);
thetaEst = zeros(noiseLen,MC,nbMethods);
aEst = zeros(noiseLen,MC,nbMethods);
bEst = zeros(noiseLen,MC,nbMethods);
timeTaken = zeros(noiseLen,nbMethods);

ii=zeros(noiseLen,MC);

CRBmat=zeros(5,noiseLen);

%Start guess for LM and fmin goes in here
[~, maxInd]=max(abs(dataNoNoise));
S0start=dataNoNoise(maxInd);
x0=startGuess(TE,TR,T1start,T2start,alpha,S0start);       
h = waitbar(0, ['Processing noiseVar 1 out of ' num2str(noiseLen) ', please wait...'],'CreateCancelBtn','disp(''Cancelled by user.'');cancelled=1;');
cancelled=0;
%% Estimation of parameters
%(can be run in parallel by parfor, MATLAB parallel toolbox required)
for k0 = 1:noiseLen
  theNoise = sqrt(noiseVar(k0)/2)*(randn(thetaLen,MC) + 1i*randn(thetaLen,MC));
  %Compute the Cramer-Rao bound
  CRBmat(:,k0)=sqrt(CRB([theta real(S0) imag(S0) a b],noiseVar(k0),deltaThetas,TE,TR));
  for modInd=1:nbMethods
    %tic;
    startTime = cputime;
    for mc = 1:MC
      %Add noise do data
      data = dataNoNoise + theNoise(:,mc);
      %Compute the estimates
      switch methodsToTest{modInd}
        case 'LORE'
          [thetaEst(k0,mc,modInd),S0Est(k0,mc,modInd),aEst(k0,mc,modInd),bEst(k0,mc,modInd)] = ...
            LORE(data,deltaThetas,TE,TR);
        case 'LORE-GN'
          [thetaInit, MInit, aInit, bInit] = ...
            LORE(data,deltaThetas,TE,TR);
          [thetaEst(k0,mc,modInd),S0Est(k0,mc,modInd),aEst(k0,mc,modInd),bEst(k0,mc,modInd),~,ii(k0,mc)] = ...
            GNORM(data,deltaThetas,[thetaInit real(MInit) imag(MInit) aInit bInit],TE,TR,1e-8,100);
        case 'LM'
          [thetaEst(k0,mc,modInd),S0Est(k0,mc,modInd),aEst(k0,mc,modInd),bEst(k0,mc,modInd)] = ...
              LMORM(data,deltaThetas,x0,TE,TR);
        case 'LMpost'
          [thetaEst(k0,mc,modInd),S0Est(k0,mc,modInd),aEst(k0,mc,modInd),bEst(k0,mc,modInd)] = ...
              LMORM(data,deltaThetas,x0,TE,TR);
              %post processing
          [S0Est(k0,mc,modInd),aEst(k0,mc,modInd),bEst(k0,mc,modInd),thetaEst(k0,mc,modInd)] = ...
            postProc(S0Est(k0,mc,modInd),aEst(k0,mc,modInd),bEst(k0,mc,modInd),thetaEst(k0,mc,modInd),TE,TR);
          otherwise
              delete(h);
              error([methodsToTest{modInd} ' is not a valid method'])
      end
    end
    timeTaken(k0,modInd) = cputime - startTime;
  end
      waitbar(k0/noiseLen, h, ['Processing noiseVar ' num2str(k0) ' out of ' num2str(noiseLen) '. Time left: ' num2str(round(sum(timeTaken(k0,:))*(noiseLen-k0))) ' s']);
    if cancelled==1
        delete(h);
        return;
    end
end
timeTakenTot = cputime - startTimeTot;
delete(h);
%% MSE calulation
thetaEstDiff = (abs((thetaEst - theta)*180/(pi)));
rmseTheta = squeeze(sqrt( mean( thetaEstDiff.^2,2 ) ));
rmseM = squeeze(sqrt( mean( abs(S0Est - S0).^2,2 ) ));
rmsea = squeeze(sqrt( mean( abs(aEst - a).^2,2 ) ));
rmseb = squeeze(sqrt( mean( abs(bEst - b).^2,2 ) ));

%Compute CRB for T1, T2 and KM0
CRB2=sqrt(CRBT1T2(T1,T2,real(KM0),imag(KM0),theta,alpha,deltaThetas,TE,TR,noiseVar));
T2Est=-TR./log(aEst);
T1Est=-TR./log((aEst*(1+ca)-bEst.*(1+aEst.^2*ca))./(aEst*(1+ca)-bEst.*(ca+aEst.^2)));
KM0Est=S0Est./(1i*sa*( 1 - exp(-TR./T1Est) )./( 1 - exp(-TR./T1Est)*ca - (exp(-TR./T1Est)-ca).*exp(-TR./T2Est).^2).*exp(-TE./T2Est));
rmseKM0 = squeeze(sqrt( mean( abs(KM0Est - KM0).^2,2 ) ));

%Remove outliers from T1 and T2 estimates for more representative plotting
outlierRateT1=zeros(noiseLen,nbMethods);
outlierRateT2=zeros(noiseLen,nbMethods);
maxDev=20; %Distance in number of standard deviations that define an outlier
rmseT1=zeros(noiseLen,nbMethods);
rmseT2=zeros(noiseLen,nbMethods);
for p1=1:noiseLen
    for p2=1:nbMethods
        T1EstTmp=T1Est(p1,:,p2);
        T2EstTmp=T2Est(p1,:,p2);
        rmseT1(p1,p2)=sqrt( mean( abs(T1EstTmp(T1EstTmp>max(0,T1-maxDev*CRB2(1,p1))&T1EstTmp<T1+maxDev*CRB2(1,p1)) - T1).^2 ) );
        rmseT2(p1,p2)=sqrt( mean( abs(T2EstTmp(T2EstTmp>max(0,T2-maxDev*CRB2(2,p1))&T2EstTmp<T2+maxDev*CRB2(2,p1)) - T2).^2 ) );
        outlierRateT1(p1,p2)=length(T1EstTmp(T1EstTmp>max(0,T1-maxDev*CRB2(1,p1))&T1EstTmp<T1+maxDev*CRB2(1,p1)))/MC;
        outlierRateT2(p1,p2)=length(T2EstTmp(T2EstTmp>max(0,T2-maxDev*CRB2(2,p1))&T2EstTmp<T2+maxDev*CRB2(2,p1)))/MC;
    end
end

%% Plotting
if plotResults
    
    figure()
    loglog(SNR_MRI, rmseT1)
    hold on;
    loglog(SNR_MRI,CRB2(1,:),'k')
    xlabel('SNR (log-scale)')
    ylabel('rMSE')
    title(['rMSE for T1 = ' num2str(T1) ', MC = ' num2str(MC)])
    legend(methodsToTest,'CRB')
    set(gca,'Xlim',[5 100],'Xtick',[5 10 15 20 30 40 60 80 100])
    fixLineStyles
    
    
    figure()
    loglog(SNR_MRI, rmseT2)
    hold on;
    loglog(SNR_MRI,CRB2(2,:),'k')
    xlabel('SNR (log-scale)')
    ylabel('rMSE')
    title(['rMSE for T2 = ' num2str(T2) ', MC = ' num2str(MC)])
    legend(methodsToTest,'CRB')
    set(gca,'Xlim',[5 100],'Xtick',[5 10 15 20 30 40 60 80 100])
    fixLineStyles
    
    figure()
    loglog(SNR_MRI, rmseKM0)
    hold on;
    loglog(SNR_MRI,CRB2(3,:),'k')
    title(['rMSE for KM0 = ' num2str(KM0) ', MC = ' num2str(MC)])
    legend(methodsToTest,'CRB')
    xlabel('SNR (log-scale)')
    ylabel('rMSE')
    set(gca,'Xlim',[5 100],'Xtick',[5 10 15 20 30 40 60 80 100])
    fixLineStyles
    
    figure()
    loglog(SNR_MRI, rmseTheta)
    hold on;
    loglog(SNR_MRI,CRBmat(1,:)*180/(pi),'k')
    title(['\theta = ' num2str(theta*180/pi) ', MC = ' num2str(MC)])
    legend(methodsToTest,'CRB')
    xlabel('SNR (log-scale)')
    ylabel('rMSE')
    set(gca,'Xlim',[5 100],'Xtick',[5 10 15 20 30 40 60 80 100])
    fixLineStyles
    
    figure()
    loglog(SNR_MRI, rmseM)
    hold on;
    loglog(SNR_MRI,sqrt(CRBmat(2,:).^2+CRBmat(3,:).^2),'k')
    title(['S_0 = ' num2str(S0) ', MC = ' num2str(MC)])
    legend(methodsToTest,'CRB')
    xlabel('SNR (log-scale)')
    ylabel('rMSE')
    set(gca,'Xlim',[5 100],'Xtick',[5 10 15 20 30 40 60 80 100])
    fixLineStyles
    
    figure()
    loglog(SNR_MRI, rmsea)
    hold on;
    loglog(SNR_MRI,CRBmat(4,:),'k')
    title(['a = ' num2str(a) ', MC = ' num2str(MC)])
    legend(methodsToTest,'CRB')
    xlabel('SNR (log-scale)')
    ylabel('rMSE')
    set(gca,'Xlim',[5 100],'Xtick',[5 10 15 20 30 40 60 80 100])
    fixLineStyles
    
    figure()
    loglog(SNR_MRI, rmseb)
    hold on;
    loglog(SNR_MRI,CRBmat(5,:),'k')
    title(['b = ' num2str(b) ', MC = ' num2str(MC)])
    legend(methodsToTest,'CRB')
    xlabel('SNR (log-scale)')
    ylabel('rMSE')
    set(gca,'Xlim',[5 100],'Xtick',[5 10 15 20 30 40 60 80 100])
    fixLineStyles

end