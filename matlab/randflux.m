function r=randflux(m,b,varargin)
% function r=randflux(m,b,varargin)

%
% $Id: randflux.m,v 1.7 2011/03/26 15:36:05 patrick Exp $
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

EPS=1e-2;

r=rand(varargin{:});

% start value is inflection point of r=CDF(x,m,b)
y = repmat(0.5*(m+sqrt(m^2+4*b^2)),size(r));

dy = (r-CDF(y,m,b))./dCDF(y,m,b);
y = y + dy;
while (max(abs(dy))>EPS),
	dy = (r-CDF(y,m,b))./dCDF(y,m,b);
	y = y + dy;
end
r=y;

function r=CDF(x,m,b)
z=(x-m)/sqrt(2)/b;
F0=F(-m/sqrt(2)/b,m,b);
Finf=m/2/b*sqrt(2*pi);
r=(F(z,m,b)-F0)/(Finf-F0);

function r=dCDF(x,m,b)
z=(x-m)/sqrt(2)/b;
F0=F(-m/sqrt(2)/b,m,b);
Finf=m/2/b*sqrt(2*pi);
r=dF(z,m,b)/(Finf-F0)/sqrt(2)/b;

function r=F(x,m,b)
r=-exp(-x.^2)+m/2/b*sqrt(2*pi)*erf(x);

function r=dF(x,m,b)
r=(2*x+m/b*sqrt(2)).*exp(-x.^2);

