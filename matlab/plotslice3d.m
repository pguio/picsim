function [y,phit,phib,phibl]=plotslice3d(filename,time,mode,seeds)
% function [y,phit,phib,phibl]=plotslice3d(filename,time,mode,seeds)

%
% $Id: plotslice3d.m,v 1.6 2011/03/26 15:36:05 patrick Exp $
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

if ~exist('seeds','var')
  seeds=[];
end

f=readhdf(filename,'potential',{time,[],[]},seeds);
fxyz=f.var;
x=f.dims{2};
y=f.dims{3};
z=f.dims{4};


g=1e-2;
Te=10;

K=g/Te;
L=max(y)-min(y);
[Y,X,Z]=ndgrid(y,x,z);
n0=K*L/(1-exp(-K*L));
f=n0.*exp(-K*Y);
y0=1/K*log(n0);
phi=-K*Te*(Y-y0);

if exist('mode','var') & mode==1,
	subplot(211)
  imagesc(x,y,squeeze(fxyz(:,:,fix(end/2))-phi(:,:,fix(end/2))))
  axis xy
  title(sprintf('%s (\\tau=%d)', tt,times));
	subplot(212)
  imagesc(y,z,squeeze(fxyz(fix(end/2),:,:)-phi(fix(end/2),:,:)))
  axis xy
  title(sprintf('%s (\\tau=%d)', f.varname, f.dims{1}));
  return
end

ps=31;
subplot(311)
size(fxyz(ps:end-ps+1,:,ps:end-ps+1))
phit=mean(mean(squeeze(fxyz(ps:end-ps+1,:,ps:end-ps+1)),3),1);
plot(y,phit) 
title(sprintf('%s (\\tau=%d)', f.varname, f.dims{1}));

bs=10;
subplot(312)
phib=mean(mean(squeeze(fxyz([1:bs end-bs+1:end],:,[1:bs end-bs+1:end])),3),1);
[p,s]=polyfit(y,phib,1);
phibl=polyval(p,y,s);
phibt=-K*Te*(y-y0);
plot(y,phib,'x',y,phibl,y,phibt)

subplot(313)
plot(y,phit-phib,'x',y,phit-phibl)
xlabel('y');
