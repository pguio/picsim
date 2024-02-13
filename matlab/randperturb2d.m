function [rx,ry]=randperturb2d(mx,sx,my,sy,alpha,varargin)
% function [rx,ry]=randperturb2d(mx,sx,my,sy,alpha,varargin)

%
% $Id: randperturb2d.m,v 1.2 2011/03/26 15:36:05 patrick Exp $
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

amax=pi/2*sx*(erf((1-mx)/sqrt(2)/sx)+erf(mx/sqrt(2)/sx))* ...
	sy*(erf((1-my)/sqrt(2)/sy)+erf(my/sqrt(2)/sy));
if (any(abs(alpha)>amax))
	error(sprintf('|alpha| must be smaller than %f',amax));
end

EPS=1e-1;

rx=rand(varargin{:});

y=mx;
%dy = (rx-CDFx(y,mx,sx,my,sy,alpha))./dCDFx(y,mx,sx,my,sy,alpha);
dy = (rx-CDFy(y,mx,sx,alpha))./dCDFy(y,mx,sx,alpha);
y = y + dy;
while (max(abs(dy))>EPS),
%	dy = (rx-CDFx(y,mx,sx,my,sy,alpha))./dCDFx(y,mx,sx,my,sy,alpha);
	dy = (rx-CDFy(y,mx,sx,alpha))./dCDFy(y,mx,sx,alpha);
	y = y + dy;
end
rx=y;

ry=rand(varargin{:});
y=my;
alpha=1/sqrt(2*pi)/sx*exp(-(rx-mx).^2/2/sx^2)*2*alpha/ ...
	(erf((1-mx)/sqrt(2)/sx)+erf(mx/sqrt(2)/sx));
dy = (ry-CDFy(y,my,sy,alpha))./dCDFy(y,my,sy,alpha);
y = y + dy;
while (max(abs(dy))>EPS),
	dy = (ry-CDFy(y,my,sy,alpha))./dCDFy(y,my,sy,alpha);
	y = y + dy;
end
ry=y;

function r=dCDFx(x,mx,sx,my,sy,alpha)
if 0
k=4*alpha/(erf((1-mx)/sqrt(2)/sx)+erf(mx/sqrt(2)/sx))/ ...
	(erf((1-my)/sqrt(2)/sy)+erf(my/sqrt(2)/sy));
beta=1+k/4*(erf((1-mx)/sqrt(2)/sx)+erf(mx/sqrt(2)/sx))* ...
	(erf((1-my)/sqrt(2)/sy)+erf(my/sqrt(2)/sy));
r=1./beta.*(1+k/(2*pi)/sx/sy.*exp(-(x-mx).^2/2/sx^2));
else
	k=2*alpha/(erf((1-mx)/sqrt(2)/sx)+erf(mx/sqrt(2)/sx));
	beta=1+k/2*(erf((1-mx)/sqrt(2)/sx)+erf(mx/sqrt(2)/sx));
	r = 1./beta.*(1+k/sqrt(2*pi)/sx.*exp(-(x-mx).^2/2/sx^2));
end

function r=CDFx(x,mx,sx,my,sy,alpha)
if 0
k=4*alpha/(erf((1-mx)/sqrt(2)/sx)+erf(mx/sqrt(2)/sx))/ ...
  (erf((1-my)/sqrt(2)/sy)+erf(my/sqrt(2)/sy));
beta=1+k/4*(erf((1-mx)/sqrt(2)/sx)+erf(mx/sqrt(2)/sx))* ...
  (erf((1-my)/sqrt(2)/sy)+erf(my/sqrt(2)/sy));
r=1./beta.*(x+k/2.*(erf((x-mx)/sqrt(2)/sx)+erf(mx/sqrt(2)/sx))/sqrt(2*pi)/sy);
else
	k=2*alpha/(erf((1-mx)/sqrt(2)/sx)+erf(mx/sqrt(2)/sx));
	beta=1+k/2*(erf((1-mx)/sqrt(2)/sx)+erf(mx/sqrt(2)/sx));
	r = 1./beta.*(x+k/2.*(erf((x-mx)/sqrt(2)/sx)+erf(mx/sqrt(2)/sx)));
end

function r=CDFy(x,m,s,alpha)
k=2*alpha/(erf((1-m)/sqrt(2)/s)+erf(m/sqrt(2)/s));
beta=1+k/2*(erf((1-m)/sqrt(2)/s)+erf(m/sqrt(2)/s));
r = 1./beta.*(x+k/2.*(erf((x-m)/sqrt(2)/s)+erf(m/sqrt(2)/s)));

function r=dCDFy(x,m,s,alpha)
k=2*alpha/(erf((1-m)/sqrt(2)/s)+erf(m/sqrt(2)/s));
beta=1+k/2*(erf((1-m)/sqrt(2)/s)+erf(m/sqrt(2)/s));
r = 1./beta.*(1+k/sqrt(2*pi)/s.*exp(-(x-m).^2/2/s^2));

