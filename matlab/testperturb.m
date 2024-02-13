function testperturb(m,s,alpha)
% function testperturb(m,s,alpha)

%
% $Id: testperturb.m,v 1.7 2011/03/26 15:36:05 patrick Exp $
%
% Copyright (c) 2000-2011 Patrick Guio <patrick.guio@gmail.com>
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

close all;

amax=sqrt(pi/2)*s*(erf((1-m)/sqrt(2)/s)+erf(m/sqrt(2)/s));
if (abs(alpha)>amax)
	error(sprintf('|alpha| must be smaller than %f',amax));
end

EPS=1e-1;

[R,CDFR] = fplot('perturbCDF',[0 1],[],[],[],m,s,alpha);

r=linspace(0.0,1.0,10);
r=CDFR;
y = m;

plot(r, y, '-', CDFR, R, '-') 
pause

dy = (r-perturbCDF(y,m,s,alpha))./dperturbCDF(y,m,s,alpha);
y = y + dy;
plot(r, y, '-', CDFR, R, '-'), pause(.1)
while (max(abs(dy))>EPS),
	dy = (r-perturbCDF(y,m,s,alpha))./dperturbCDF(y,m,s,alpha);
	y = y + dy;
	plot(r, y, '-', CDFR, R, '-')
%	plot(r,abs(dy./y))
	pause(.1)
end

err=abs(r-perturbCDF(y,m,s,alpha));
fprintf(1,'mean(error)=%.2e std(error)=%.2e\n', mean(err), std(err));

ii=find(err~=0);
subplot(211), plot(r, y, '-', CDFR, R, '-')
subplot(212), plot(log(err(ii)), y(ii), '-')

function r=perturbCDF(x,m,s,alpha)
k = 2*alpha/(erf((1-m)/sqrt(2)/s) + erf(m/sqrt(2)/s));
beta = 1 + k/2*(erf((1-m)/sqrt(2)/s) + erf(m/sqrt(2)/s));
r = 1/beta*(x + k/2*(erf((x-m)/sqrt(2)/s) + erf(m/sqrt(2)/s)));

function r=dperturbCDF(x,m,s,alpha)
k = 2*alpha/(erf((1-m)/sqrt(2)/s) + erf(m/sqrt(2)/s));
beta = 1 + k/2*(erf((1-m)/sqrt(2)/s) + erf(m/sqrt(2)/s));
r = 1/beta*(1 + k/sqrt(2*pi)/s*exp(-(x-m).^2/2/s^2));

function r=dF2(x,m,s,alpha)
k = 2*alpha/(erf((1-m)/sqrt(2)/s) + erf(m/sqrt(2)/s));
beta = 1 + k/2*(erf((1-m)/sqrt(2)/s) + erf(m/sqrt(2)/s));
r = 1/beta*(-k/sqrt(2*pi)/s^3*(x-m)*exp(-(x-m).^2/2/s^2));

