function tracer=tracer_init(particles,scheduler)
% function tracer=tracer_init(particles,scheduler)

%
% $Id: tracer_init.m,v 1.5 2011/03/26 15:36:06 patrick Exp $
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

ii = find(sum(particles.R.^2) < 4 & abs(particles.R(1,:)) < 0.3 & ...
          sum(particles.V.^2) < 0.05 & ...
          abs(sum(particles.V([2 3],:).^2)./particles.V(1,:)) < 0.2 & ...
          abs(sum(particles.V([2 3],:).^2)./particles.V(1,:)) < 0.3)
tracer.index = ii(1:2);
tracer.X=zeros(scheduler.maxiter,length(tracer.index));
tracer.Y=zeros(scheduler.maxiter,length(tracer.index));
if scheduler.dimr==3,
	tracer.Z=zeros(scheduler.maxiter,length(tracer.index));
end
