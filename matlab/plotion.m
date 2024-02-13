function [is,ks,fs]=plotion(filename, subdim, n_av, plotcmd, seeds)
% function [is,ks,fs]=plotion(filename, subdim, n_av, plotcmd, seeds)

%
% $Id: plotion.m,v 1.8 2011/03/26 15:36:05 patrick Exp $
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

f=readhdf(filename,'probes',subdim,seeds);
t=f.dims{1};
x=f.dims{2};
ftx=f.var;

if 1
	if 0
		t=t(1:2:end-1);
		ftx=ftx(:,1:2:end-1);
	else
		t=t(2:2:end);
		ftx=ftx(:,2:2:end);
	end
end

if rem(size(ftx,2), n_av)~=0
	error(sprintf('n_av must be a multiple of %d',size(ftx,2)))
end

ftk=fft(ftx);
ks=fftshift(2*pi*[-size(ftx,1)/2:size(ftx,1)/2-1]'/(x(end)-x(1)));

ftk=ftk(1:fix(end/2),:);
ks=ks(1:fix(end/2));

is=zeros(size(ftk,1), size(ftk,2)/n_av);
for ik=1:size(ftk,1),
	it=[1:size(ftk,2)/n_av];
	specs=zeros(n_av,size(ftk,2)/n_av);
	for k=1:n_av,
		% fprintf(1,'min(it)=%d max(it)=%d\n',min(it),max(it));
		specs(k,:) = fftshift(abs(fft(conj(ftk(ik,it))')).^2)';
		it = it + size(ftk,2)/n_av;
	end
	is(ik,:) = mean(specs);
end
fs=1/(t(2)-t(1))*[-size(is,2)/2:size(is,2)/2-1]/size(is,2);

%return

iks=2:8:length(ks);
iks=3:3:42;

n=length(ks(iks));
n2=fix(length(ks(iks))/2);

if ~exist('plotcmd')
  plotcmd='plot';
end

te=10;
ti=1;
method=3;
switch method,
	case 1
		cs = sqrt(te+ti);
		fis = ks*cs/(2*pi);
	case 2
		cs = 1/sqrt(2)*sqrt(1+sqrt(1+12*ti/te))*sqrt(te);
		fis = ks*cs/(2*pi);
	case 3
		omi = 1;
		ke = sqrt(0.1);
		fis = ks*omi/sqrt(2)/ke./sqrt(1+ks.^2/ke^2).* ...
			sqrt(1+12*ti/te*sqrt(1+ks.^2/ke^2))/(2*pi);
end

for i=1:n2,
	subplot(n2,2,2*i-1)
	if ks(i)<0.0,
		feval(plotcmd,-fs,is(iks(i),:))
	else
		feval(plotcmd,fs,is(iks(i),:))
	end
	h=line([fis(iks(i)) fis(iks(i))],get(gca,'ylim')); 
	set(h,'linestyle','--');
	h=line(-[fis(iks(i)) fis(iks(i))],get(gca,'ylim'));
	set(h,'linestyle','--');
	title(sprintf('<|n(k,f)|^2> (k=%.2f)',ks(iks(i))));

	subplot(n2,2,2*i)
	if ks(iks(i+n2))<0.0,
		feval(plotcmd,-fs,is(iks(i+n2),:))
	else
		feval(plotcmd,fs,is(iks(i+n2),:))
	end
	h=line([fis(iks(i+n2)) fis(iks(i+n2))],get(gca,'ylim')); 
	set(h,'linestyle','--');
	h=line(-[fis(iks(i+n2)) fis(iks(i+n2))],get(gca,'ylim'));
	set(h,'linestyle','--');
	title(sprintf('<|n(k,f)|^2> (k=%.2f)',ks(iks(i+n2))));
end
