function y=fluxCDF(x,m,b)
% function y=fluxCDF(x,m,b)

%
% $Id: fluxCDF.m,v 1.3 2011/03/26 15:36:05 patrick Exp $
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

z = (x-m)/sqrt(2)/b;
F0 = F(-m/sqrt(2)/b,m,b);
Finf = m/2/b*sqrt(2*pi);
y = (F(z,m,b)-F0)/(Finf-F0);

function r=F(x,m,b)
r = -exp(-x.^2) + m/b/2*sqrt(2*pi)*erf(x);

