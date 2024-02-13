function testdriven2d
% function testdriven2d

%
% $Id: testdriven2d.m,v 1.7 2011/03/26 15:36:05 patrick Exp $
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

m=0.5;
s=0.1;
alpha=0.2;

maxiter=200;
nbins=20;

Npart=50000;
dt=0.2;
omega=4*pi/(maxiter-1)/dt;

x=linspace(0,1,nbins);

fxt=zeros(maxiter,length(x));
time=linspace(0,dt*(maxiter-1),maxiter);
for i=1:maxiter
	a = -alpha/(omega*dt)*(cos(omega*(time(i)+dt))-cos(omega*time(i)));
	r=randperturb(m, s, a, 1, fix(Npart*(1+a)));
	if 1
		fxt(i,:)=histc(r,x);
	else
		ix=fix(r*(nbins-1))+1;
		for k=1:length(ix),
			fxt(i,ix(k))=fxt(i,ix(k))+1;
		end
	end
	fprintf(1,'.');
end
fprintf(1,'\n');

x=x(1:end-1);
fxt=fxt(:,1:end-1);

subplot(211), 
imagesc(time,x,fxt')
axis xy
	
x=repmat(x,maxiter,1);
ft=d01gaf(x',fxt');

subplot(212)
plot(time, ft)
fprintf(1,'amplitude=%f (alpha=%f)\n',sqrt(2)*std(ft)/mean(ft),alpha);

