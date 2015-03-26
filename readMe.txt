Short description of the contents of the LORE-GN package. The code comes with no warranty.

Please refer to the paper for more details: 
"Parameter estimation approach to banding artifact reduction in balanced steady-state free precession"
http://onlinelibrary.wiley.com/doi/10.1002/mrm.24986/full

Written by: Marcus Bjork, R. Reeve Ingle, Erik Gudmundson, Joelle K. Barral, 2013 
Contact: marcus.bjork.85@gmail.com

Copyright, Marcus Bjork, Uppsala University, 2013

This file is part of the LORE-GN package.

LORE-GN is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

LORE-GN is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with LORE-GN.  If not, see <http://www.gnu.org/licenses/>.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Runnable scripts:
runLOREGN (Figs. 4, 5, 6, 7 and 8)
MSEvsSNR (Figs. 1 and 2)
CRBplot (Fig. 3)

Main functions:
LORE     	- LORE algorithm
GNORM    	- Gauss-Newton implementation
postProc 	- Post processing step

Help functions:
SNRcalc		- Estimates SNR
mat2image	- Makes a big image from 3D matrix of images (for plotting)
growShrink	- Shrink or grow a mask
genData		- Generate data based on the bSSFP model
genDataMat	- Same as genData but with extended matrix support (slower)
CRB		- Computes CRB based on S_0, a, b, theta parameterization
CRBT1T2		- Computes CRB based on KM_0, T1, T2, theta parameterization
fixLineStyles	- Sets the line styles for the MSE plots to correspond with the paper

For comparison:
LMORM		- Uses Levenberg-Marquardt to estimate parameters
startGuess	- Computes starting guess vector for LMORM 
SoS		- Sum-of-Squares method
maxIntensity	- Maximum Intensity method
XSORM		- Cross-Solution from Xiang et.al.