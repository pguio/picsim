function [y,phit,phib,phibl]=plotslice2d(filename,time,mode,seeds)
% function [y,phit,phib,phibl]=plotslice2d(filename,time,mode,seeds)

%
% $Id: plotslice2d.m,v 1.8 2011/03/26 15:36:05 patrick Exp $
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
fxy=f.var;
x=f.dims{2};
y=f.dims{3};

g=1e-2;
Te=10;

K=g/Te;
L=max(y)-min(y);
[Y,X]=ndgrid(y,x);
n0=K*L/(1-exp(-K*L));
f0=n0.*exp(-K*Y);
y0=1/K*log(n0);
phi=-K*Te*(Y-y0);

if exist('mode','var') & mode==1,
	subplot(111)
	imagesc(x,y,fxy-phi)
	axis xy
	title(sprintf('%s (\\tau=%d)', f.varname, f.dims{1}));
	colorbar
	return
end

ps=40;
subplot(311)
phit=mean(fxy(:,ps:end-ps+1)');
plot(y,phit) 
title(sprintf('%s (\\tau=%d)', f.varname, f.dims{1}));

bs=5;
subplot(312)
phib=mean(fxy(:,[1:bs end-bs+1:end])');
[p,s]=polyfit(y,phib,1);
phibl=polyval(p,y,s);
phibt=-K*Te*(y-y0);
plot(y,phib,'x',y,phibl,y,phibt)

subplot(313)
plot(y,phit-phib,'x',y,phit-phibl,y,phit-phibt)
xlabel('y');
