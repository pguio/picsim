function checkperturb(m,s,alpha)
% function checkperturb(m,s,alpha)

%
% $Id: checkperturb.m,v 1.9 2011/03/26 15:36:05 patrick Exp $
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

x=linspace(0,1,100);

[y,y0]=perturb(x,m,s,alpha);
plot(x,y,x,y0);

i=d01gaf(x(:),y(:));
i0=d01gaf(x(:),y0(:));

fprintf('i=%f\n', i);
fprintf('i0=%f\n', i0);
fprintf('(i-i0)/i0=%f\n', (i-i0)/i0);

function [y,y0]=perturb(x,m,s,alpha)

amax=sqrt(pi/2)*s*(erf((1-m)/sqrt(2)/s)+erf(m/sqrt(2)/s));
if (abs(alpha)>amax)
	error(sprintf('|alpha| must be smaller than %f',amax));
end
k=2*alpha/(erf((1-m)/sqrt(2)/s)+erf(m/sqrt(2)/s));
beta=1+k/2*(erf((1-m)/sqrt(2)/s)+erf(m/sqrt(2)/s));
y = 1/beta*(1 + k/sqrt(2*pi)/s*exp(-(x-m).^2/2/s^2));
y0 = repmat(1/beta,size(y));

if 0

if alpha>=0
	k=2*alpha/(erf((1-m)/sqrt(2)/s)+erf(m/sqrt(2)/s));
	beta=1+k/2*(erf((1-m)/sqrt(2)/s)+erf(m/sqrt(2)/s));
	y = 1/beta*(1 + k/sqrt(2*pi)/s*exp(-(x-m).^2/2/s^2));
	y0 = repmat(1/beta,size(y));
else
	k=2*alpha/(-2*alpha/sqrt(2*pi)/s+erf((1-m)/sqrt(2)/s)+erf(m/sqrt(2)/s));
	beta=1-k/sqrt(2*pi)/s+k/2*(erf((1-m)/sqrt(2)/s)+erf(m/sqrt(2)/s));
	y = 1/beta*(1 - k/sqrt(2*pi)/s*(1-exp(-(x-m).^2/2/s^2)));
	y0 = repmat(1/beta*(1 - k/sqrt(2*pi)/s),size(y));
	(beta-(1 - k/sqrt(2*pi)/s))/(1+k/sqrt(2*pi)/s)
end

end
