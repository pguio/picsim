function r=randperturb(m,s,alpha,varargin)
% function r=randperturb(m,s,alpha,varargin)

%
% $Id: randperturb.m,v 1.8 2011/03/26 15:36:05 patrick Exp $
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

amax=sqrt(pi/2)*s*(erf((1-m)/sqrt(2)/s)+erf(m/sqrt(2)/s));
if (any(abs(alpha)>amax))
	error(sprintf('|alpha| must be smaller than %f',amax));
end

EPS=1e-1;

r=rand(varargin{:});

y=m;
dy = (r-CDF(y,m,s,alpha))./dCDF(y,m,s,alpha);
y = y + dy;
while (max(abs(dy))>EPS),
	dy = (r-CDF(y,m,s,alpha))./dCDF(y,m,s,alpha);
	y = y + dy;
end
r=y;

function r=dCDF(x,m,s,alpha)
r = dF(x,m,s,alpha);

function r=CDF(x,m,s,alpha)
r = F(x,m,s,alpha);

function r=F(x,m,s,alpha)
k=2*alpha/(erf((1-m)/sqrt(2)/s)+erf(m/sqrt(2)/s));
beta=1+k/2*(erf((1-m)/sqrt(2)/s)+erf(m/sqrt(2)/s));
r = 1./beta.*(x+k/2.*(erf((x-m)/sqrt(2)/s)+erf(m/sqrt(2)/s)));

function r=dF(x,m,s,alpha)
k=2*alpha/(erf((1-m)/sqrt(2)/s)+erf(m/sqrt(2)/s));
beta=1+k/2*(erf((1-m)/sqrt(2)/s)+erf(m/sqrt(2)/s));
r = 1./beta.*(1+k/sqrt(2*pi)/s.*exp(-(x-m).^2/2/s^2));

