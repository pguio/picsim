function [kx,ky,dkx,dky]=findkxky(x,y,z)
% function [kx,ky,dkx,dky]=findkxky(x,y,z)

%
% $Id: findkxky.m,v 1.2 2011/03/26 15:36:05 patrick Exp $
%
% Copyright (c) 2001-2011 Patrick Guio <patrick.guio@gmail.com>
% All Rights Reserved.
%
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2.  of the License, or (at your
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.

nx=length(x);
ny=length(y);

lx=max(x)-min(x);
ly=max(y)-min(y);

nxfft=256;
nyfft=512;

kxs=2*pi*[-nxfft/2:nxfft/2-1]'/(nxfft-1)*(nx-1)/lx;
kys=2*pi*[-nyfft/2:nyfft/2-1]'/(nyfft-1)*(ny-1)/ly;

subplot(121)
imagesc(x,y,z);

subplot(122)
ftz=fftshift(abs(fft2(z-mean(z(:)),nyfft,nxfft)).^2);
imagesc(kxs,kys,log(ftz))

mxz=max(ftz(:));
[i,j]=find(ftz==mxz);

fprintf(1,'Calculated values for kx and ky\n');
fprintf(1,'kx=%+.6f, ky=%+.6f\n', [kxs(j), kys(i)]');
fprintf(1,'dkx=%f, dky=%f\n', max(diff(kxs)), max(diff(kys)))

kx=kxs(j);
ky=kys(i);

dkx= max(diff(kxs));
dky=max(diff(kys));
