function r=perturbCDF(x,m,s,alpha)
% function r=perturbCDF(x,m,s,alpha)

%
% $Id: perturbCDF.m,v 1.4 2011/03/26 15:36:05 patrick Exp $
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
if (abs(alpha)>amax)
	error(sprintf('|alpha| must be smaller than %f',amax));
end
k=2*alpha/(erf((1-m)/sqrt(2)/s)+erf(m/sqrt(2)/s));
beta=1+k/2*(erf((1-m)/sqrt(2)/s)+erf(m/sqrt(2)/s));
r = 1/beta*(x+k/2*(erf((x-m)/sqrt(2)/s)+erf(m/sqrt(2)/s)));

