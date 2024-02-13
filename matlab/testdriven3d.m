function testdriven3d
% function testdriven3d

%
% $Id: testdriven3d.m,v 1.12 2011/03/26 15:36:05 patrick Exp $
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

mx=0.5;
sx=0.1;
my=0.5;
sy=0.1;
%alpha=0.0995;
alpha=0.06;

maxiter=200;
nbins=20;

Npart=100000;
dt=0.2;
omega=4*pi/(maxiter-1)/dt;

x=linspace(0,1,nbins);
y=linspace(0,1,nbins);

fxt=zeros(maxiter,length(x),length(y));
ft=zeros(maxiter,1);
time=linspace(0,dt*(maxiter-1),maxiter);
for i=1:maxiter
	if maxiter==200,
		a = -alpha/(omega*dt)*(cos(omega*(time(i)+dt))-cos(omega*time(i)));
	else
		a = alpha;
	end
	if 1
	rx=randperturb(mx, sx, a, 1, fix(Npart*(1+a)));
	ay=1/sqrt(2*pi)/sx*exp(-(rx-mx).^2/2/sx^2)*2/ ...
		(erf((1-mx)/sqrt(2)/sx)+erf(mx/sqrt(2)/sx));
%	fprintf(1,'a=%f ay=[%f, %f]\n', a, min(ay), max(ay));
	ry=randperturb(my, sy, ay*a, 1, fix(Npart*(1+a)));
	else
		[rx,ry]=randperturb2d(mx, sx, my, sy, a, 1, fix(Npart*(1+a)));
	end
	ix=fix(rx*(nbins-1))+1;
	jy=fix(ry*(nbins-1))+1;
	for k=1:length(ix),
		fxt(i,ix(k),jy(k))=fxt(i,ix(k),jy(k))+1;
	end
	ft(i)=d01gaf2d(x(1:end-1)',y(1:end-1),squeeze(fxt(i,1:end-1,1:end-1))');
	fprintf(1,'.');
	subplot(111), 
	imagesc(x(1:end-1),y(1:end-1),squeeze(fxt(i,1:end-1,1:end-1)));
	axis xy
	drawnow
end
fprintf(1,'\n');

if  maxiter==200
	subplot(111)
	plot(time, ft)
	fprintf(1,'amplitude=%f (alpha=%f)\n',sqrt(2)*std(ft)/mean(ft),alpha);
else
	subplot(111)
	imagesc(x(1:end-1),y(1:end-1),squeeze(mean(fxt(:,1:end-1,1:end-1),1)))
end
