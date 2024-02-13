function [is,ks,fs]=plotiongravity(filename, n_av, plotcmd, seeds)
% function [is,ks,fs]=plotiongravity(filename, n_av, plotcmd, seeds)

%
% $Id: plotiongravity.m,v 1.6 2011/03/26 15:36:05 patrick Exp $
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

f=readhdf(filename,'probes',[],seeds);
t=f.dims{1};
x=f.dims{2};
ftx=f.var;

if rem(size(ftx,2), n_av)~=0
	error(sprintf('n_av must be a multiple of %d',size(ftx,2)))
end

g=1e-2;
Te=10;

K=g/Te;
L=max(x)-min(x);
n0=K*L/(1-exp(-K*L));
f=n0.*exp(-K*x);
x0=1/K*log(n0);
phi=repmat(-K*Te*(x'-x0),1,length(t));

if 1
	subplot(211), imagesc(t,x,ftx), axis xy, colorbar
	subplot(212), imagesc(t,x,ftx-phi), axis xy, colorbar
	pause
end

ftx=ftx-phi;



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

iks=2:4:length(ks);

n=length(ks(iks));
n2=fix(length(ks(iks))/2);

if ~exist('plotcmd')
  plotcmd='plot';
end

for i=1:n2,
	subplot(n2,2,2*i-1)
	if ks(i)<0.0,
		feval(plotcmd,-fs,is(iks(i),:))
	else
		feval(plotcmd,fs,is(iks(i),:))
	end
	h=line([ks(iks(i)) ks(iks(i))]/(2*pi),get(gca,'ylim')); 
	set(h,'linestyle','--');
	h=line(-[ks(iks(i)) ks(iks(i))]/(2*pi),get(gca,'ylim'));
	set(h,'linestyle','--');
	title(sprintf('<|n(k,f)|^2> (k=%.2f)',ks(iks(i))));

	subplot(n2,2,2*i)
	if ks(iks(i+n2))<0.0,
		feval(plotcmd,-fs,is(iks(i+n2),:))
	else
		feval(plotcmd,fs,is(iks(i+n2),:))
	end
	h=line([ks(iks(i+n2)) ks(iks(i+n2))]/(2*pi),get(gca,'ylim')); 
	set(h,'linestyle','--');
	h=line(-[ks(iks(i+n2)) ks(iks(i+n2))]/(2*pi),get(gca,'ylim'));
	set(h,'linestyle','--');
	title(sprintf('<|n(k,f)|^2> (k=%.2f)',ks(iks(i+n2))));
end
