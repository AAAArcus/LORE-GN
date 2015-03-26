%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate off-resonance and corrected (band supressed) image using 3 or
% more phase cycled complex-valued images obtained with bSSFP using
% arbitrary TE and TR.
%
% Data must include a 3D data set (im) with phases along the third dimension
% and a vector with the corresponding phases (deltaThetas). To estimate T1 and
% KM0 the flip angle (alpha) is needed. It is possible to use an image of
% estimated alphas of the same size as the images used as data.
%
% im - 3D double matrix. Example: im = im(y_spaceCoord, x_spaceCoord, deltaThetas)
% deltaThetas - phase increment vector in radians. Example: deltaThetas = [0 pi/2 pi 3*pi/2]'
% alpha - flip angle in degrees. Either scalar (assumed) or matrix (measured).
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

%% User input (uncomment selected dataset)
%dataSet='phantom_1_5T'; %Figure 4
dataSet = 'brain_1_5T'; %Figure 5, 6, 7
%dataSet='brain_7T'; % Figure 8

%Uncomment method of choice
%setAlgorithm = 'LORE';      %LORE : linearization, biased estimate used for initialization
setAlgorithm = 'LORE-GN';   %LORE-GN :  linearization then custom Gauss-Newton algorithm
%setAlgorithm = 'LM';         %LM : Levenberg-Marquardt, lsqnonlin MATLAB function with fixed startguess
%setAlgorithm = 'LMpost';     %LM : Levenberg-Marquardt, lsqnonlin MATLAB function with fixed startguess and the proposed post processing

plotData=1;
if strcmp(dataSet(end-3:end),'1_5T')
    maskFactor = 0.15;  % procentage treshold for noise
else
    maskFactor = 0.06;
end
maxIter=100;    %Maximum number of iteration in GN
tol=1e-8;  %Maximum norm of gradient for GN minimum


%% Prel
loadPath = '../Data/';
load([loadPath dataSet])

deltaThetas=deltaThetas(:);

%Get image size
[nbrow, nbcol, nbrImages] = size(im);

%Normalize for numerical stability
imageSoS=SoS(im);
im=im/max(imageSoS(:));
imageSoS=imageSoS/max(imageSoS(:));

%% Mask the background and estimate SNR 
%(points dimmer than maskFactor*[the brightest point]) are removed
mask=(imageSoS > maskFactor*max( imageSoS(:) ));

%Plot the mask
if plotData
  figure,
  imshow(abs(mask(:,:)),[])
  title(['Mask for data set ' dataSet],'Interpreter','none')
  pause(0.01) %To make sure the mask figure is displayed prior to calculations
end

%Index for non-masked voxels in image
maskInds = find(mask);

%SNR
[~, SNR_MRI, SNRim, imNoise]=SNRcalc(im,mask);
display(['SNR: ' num2str(round(SNR_MRI))])


%% Find the off-resonance and image
%Number of voxels to process
nVoxAll = length(maskInds);

%Pre-allocation
thetaEstTmp = zeros(nVoxAll, 1);
S0EstTmp = zeros(nVoxAll, 1);
aEstTmp = zeros(nVoxAll, 1);
bEstTmp = zeros(nVoxAll, 1);
lossTmp = zeros(nVoxAll, 1);
nbrIter=zeros(nVoxAll, 1);

%Start guess for LM
x0=startGuess(TE,TR,500,100,alpha,1);

%Print info in MATLAB command window
disp(['Method: ' setAlgorithm])
%disp(['Processing ' num2str(nVoxAll) ' voxels.']);
h = waitbar(0, sprintf('Processing %d voxels', nVoxAll),'CreateCancelBtn','disp(''Cancelled by user.'');cancelled=1;');

if strcmp(dataSet,'brain_7T') && strcmp(setAlgorithm,'LORE-GN')
    button=questdlg(sprintf('%s \n%s','The 7T dataset has too low SNR for stable execution of Gauss-Newton (matrix inversion). It is recommeded to use LORE alone in this case.','Do you want do continue anyway (results may be inaccurate)?'),'Warning','Yes','No','Yes');
    if strcmp(button,'No')
        disp('Execution aborted by user')
        return;
    end
end
cancelled=0;

tic;
%Loop over all pixels (can be run in parallel by parfor, MATLAB parallel toolbox required)
for jj = 1:nVoxAll
   
    [rowind, colind]=ind2sub([nbrow nbcol],maskInds(jj));
    imagePoint = squeeze(im(rowind,colind,:));  %Measurements
    switch(setAlgorithm)
      case 'LORE'
        [thetaEstTmp(jj),S0EstTmp(jj),aEstTmp(jj),bEstTmp(jj)] = ...
          LORE(imagePoint,deltaThetas,TE,TR);
      case 'LORE-GN'
        [thetaEstTmp(jj),S0EstTmp(jj),aEstTmp(jj),bEstTmp(jj)] = ...
          LORE(imagePoint,deltaThetas,TE,TR);
        [thetaEstTmp(jj),S0EstTmp(jj),aEstTmp(jj),bEstTmp(jj),lossVec, nbrIter(jj)] = ...
          GNORM(imagePoint,deltaThetas,[thetaEstTmp(jj),real(...
          S0EstTmp(jj)),imag(S0EstTmp(jj)),real(aEstTmp(jj)),bEstTmp(jj)],TE,TR,tol,maxIter);
      case 'LM'
        [thetaEstTmp(jj),S0EstTmp(jj),aEstTmp(jj),bEstTmp(jj)] = ...
          LMORM(imagePoint,deltaThetas,x0,TE,TR);
      case 'LMpost'
        [thetaEstTmp(jj),S0EstTmp(jj),aEstTmp(jj),bEstTmp(jj)] = ...
          LMORM(imagePoint,deltaThetas,x0,TE,TR);
        [S0EstTmp(jj),aEstTmp(jj),bEstTmp(jj),thetaEstTmp(jj)]=postProc(S0EstTmp(jj),aEstTmp(jj),bEstTmp(jj),thetaEstTmp(jj),TE,TR);
      otherwise
        delete(h);
        error(['Alg. ' setAlgorithm ' is unknown!'])
    end
    %Compute criterion
    dataEst=genData(S0EstTmp(jj),aEstTmp(jj),bEstTmp(jj),thetaEstTmp(jj),deltaThetas,TE,TR);
    lossTmp(jj)=1/norm(imagePoint)*sum(abs(imagePoint-dataEst).^2);
    
    %Update waitbar every 100 pixels
    if rem(jj,100)==0 || jj==nVoxAll
        timeLeft=toc/jj*(nVoxAll-jj);
        waitbar(jj/nVoxAll, h, sprintf('Processing %d voxels, %g percent done, time left: %d s',nVoxAll,round(100*jj/nVoxAll),round(timeLeft)));
    end
    if cancelled==1
        delete(h);
        return;
    end
end
timeTaken = toc;
delete(h);

%Go back to nbrow x nbcol image parameter estimates
thetaEst = zeros(nbrow,nbcol);
thetaEst(maskInds) = thetaEstTmp;
S0Est = zeros(nbrow,nbcol);
S0Est(maskInds) = S0EstTmp;
aEst = zeros(nbrow,nbcol);
aEst(maskInds) = aEstTmp;
bEst = zeros(nbrow,nbcol);
bEst(maskInds) = bEstTmp;

%Compute T1 T2 and KM0 estimates (suboptimal)
ca=cos(alpha*pi/180);
sa=cos(alpha*pi/180);
T2Est = zeros(nbrow,nbcol);
T2EstTmp=-TR./log(aEstTmp);
T2Est(maskInds)=real(T2EstTmp); 
T1Est = zeros(nbrow,nbcol);
T1EstTmp=-TR./log((aEstTmp.*(1+ca)-bEstTmp.*(1+aEstTmp.^2.*ca))./(aEstTmp.*(1+ca)-bEstTmp.*(ca+aEstTmp.^2)));
T1Est(maskInds)=real(T1EstTmp);
KM0Est=S0Est./(1i*sa*( 1 - exp(-TR./T1Est) )./( 1 - exp(-TR./T1Est)*ca - (exp(-TR./T1Est)-ca).*exp(-TR./T2Est).^2).*exp(-TE./T2Est));

%Loss image
loss=inf*ones(nbrow,nbcol);
loss(maskInds)=lossTmp;

%LORE-GN reconstruction(s)
reconThetas=[0 pi/4 pi/2 3*pi/4 pi];
imRecons=abs(genDataMat(S0Est,aEst,bEst,0,reconThetas,TE,TR));
imRecon=imRecons(:,:,3);

disp(['Processed ' num2str(nVoxAll) ' voxels in ' num2str(round(timeTaken)) ' sec (' num2str(timeTaken/60) ' mins)']);

%Comparing reconstructions
imageMI=maxIntensity(im);
imageXS=abs(XSORM(im));

%Mask for plotting XS MI and SoS
maskNAN=+mask;
maskNAN(mask==0)=nan;

%% Display results
if plotData
    
    %Get a suitable thershold for display
    tmp=sort(reshape(abs(mat2image(im)),[],1));
    maxLim=tmp(end-round(0.02*nVoxAll));
        
    figure, imshow(abs(mat2image(im)),[0 maxLim])
    title('Original data (magnitude)')
    set(gcf,'position',[500 500 4*256 285])
    set(gca,'position',[0 0 1 256/285])
    
    figure, imshow(thetaEst*1000/(2*pi*TR),[-100 100])
    h=colorbar;
    ylabel(h,'[Hz]')
    set(gcf,'position',[500 500 360 256/0.95])
    set(gca,'position',[0 0.025 256/340 0.95])
    
    if exist('B0map','var')
        figure,
        imshow(B0map,[-150 150])
        title('Measured \theta-map (Hz)')
        h=colorbar;
        ylabel(h,'[Hz]')
        set(gcf,'position',[500 500 360 256/0.95])
        set(gca,'position',[0 0.025 256/340 0.95])
    end
    
    figure, imshow(imRecon,[0 maxLim])
    title('Reconstruction with \theta=\pi/2')
    set(gcf,'position',[500 500 256 285])
    set(gca,'position',[0 0 1 256/285])
    if ~strcmp(dataSet,'brain_7T')
        figure, imshow(mat2image(imRecons),[0 maxLim])
        title('Reconstructions at different \theta')
        set(gcf,'position',[500 500 length(reconThetas)*256 285])
        set(gca,'position',[0 0 1 256/285])
    end
    
    
    figure, imshow(maskNAN.*imageSoS,[0 maxLim])
    title('Sum of Squares (SoS)')
    set(gcf,'position',[500 500 256 285])
    set(gca,'position',[0 0 1 256/285])
   
    figure, imshow(maskNAN.*imageMI,[0 maxLim])
    title('Maximum Intensity (MI)')
    set(gcf,'position',[500 500 256 285])
    set(gca,'position',[0 0 1 256/285])
    
    figure, imshow(maskNAN.*imageXS,[0 maxLim])
    title('Cross solution (XS)')
    set(gcf,'position',[500 500 256 285])
    set(gca,'position',[0 0 1 256/285])
    
    %Plot the difference images
    if strcmp(dataSet,'brain_1_5T') && strcmp(setAlgorithm(1:4),'LORE')
        a=imageSoS(mask)'*imRecon(mask)/(imageSoS(mask)'*imageSoS(mask));
        figure,imshow(-abs(imRecon)+a*maskNAN.*imageSoS,[-0.02 0.02])
        title('Difference with SoS')
        set(gcf,'position',[500 500 256 285])
        set(gca,'position',[0 0 1 256/285])

        a=imageMI(mask)'*imRecon(mask)/(imageMI(mask)'*imageMI(mask));
        figure,imshow(-abs(imRecon)+a*maskNAN.*imageMI,[-0.02 0.02])
        title('Difference with MI')
        set(gcf,'position',[500 500 256 285])
        set(gca,'position',[0 0 1 256/285])

        a=imageXS(mask)'*imRecon(mask)/(imageXS(mask)'*imageXS(mask));
        figure,imshow(-abs(imRecon)+a*maskNAN.*imageXS,[-0.02 0.02])
        title('Difference with XS')
        set(gcf,'position',[500 500 256 285])
        set(gca,'position',[0 0 1 256/285])
    end
end


