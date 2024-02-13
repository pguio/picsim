function testflux(m,b)
% function testflux(m,b)

%
% $Id: testflux.m,v 1.8 2011/03/26 15:36:05 patrick Exp $
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

EPS=1e-1;

[R,CDFR] = fplot('fluxCDF',[0 40],[],[],[],m,b);

r=linspace(0.0,1.0,10);
r=CDFR;

% start value is inflection point of fluxCDF(y,m,b)
y = repmat(0.5*(m+sqrt(m^2+4*b^2)),size(r));

plot(r, y, '-', CDFR, R, '-') 
pause;

dy = (r-fluxCDF(y,m,b))./dfluxCDF(y,m,b);
y = y + dy;
plot(r, y, '-', CDFR, R, '-'), pause(.1)
while (max(abs(dy))>EPS),
	dy = (r-fluxCDF(y,m,b))./dfluxCDF(y,m,b);
	y = y + dy;
	plot(r, y, '-', CDFR, R, '-')
%	plot(r,abs(dy./y))
	pause(.1)
end

err=abs(r-fluxCDF(y,m,b));
fprintf(1,'mean(error)=%.2e std(error)=%.2e\n', mean(err), std(err));

ii=find(err~=0);
subplot(211), plot(r, y, '-', CDFR, R, '-')
subplot(212), plot(log(err(ii)), y(ii), '-')

function r=fluxCDF(x,m,b)
z = (x-m)/sqrt(2)/b;
F0 = F(-m/sqrt(2)/b,m,b);
Finf = m/2/b*sqrt(2*pi);
r = (F(z,m,b)-F0)/(Finf-F0);

function r=dfluxCDF(x,m,b)
z = (x-m)/sqrt(2)/b;
F0 = F(-m/sqrt(2)/b,m,b);
Finf = m/2/b*sqrt(2*pi);
r = dF(z,m,b)/(Finf-F0)/sqrt(2)/b;

function r=F(x,m,b)
r = -exp(-x.^2) + m/b/2*sqrt(2*pi)*erf(x);

function r=dF(x,m,b)
r = (2*x + m/b*sqrt(2)).*exp(-x.^2);

function r=d2fluxCDF(x,m,b)
z = (x-m)/sqrt(2)/b;
F0 = F(-m/sqrt(2)/b,m,b);
Finf = m/2/b*sqrt(2*pi);
r = d2F(z,m,b)/(Finf-F0)/2/b^2;

function r=d2F(x,m,b)
r = -2*exp(-x.^2).*(-b+2*x.^2*b+x*m*sqrt(2))/b;


