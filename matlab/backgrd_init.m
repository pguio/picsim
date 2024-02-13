function particles=backgrd_init(particles,scheduler)
% function particles=backgrd_init(particles,scheduler)

%
% $Id: backgrd_init.m,v 1.5 2011/03/26 15:36:05 patrick Exp $
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

% Thermal temperature
particles.vth=sqrt(particles.temp/particles.mass);

% Limits of simulation domain
particles.rmin=scheduler.rmin;
particles.rmax=scheduler.rmax;
particles.deltar=particles.rmax-particles.rmin;

% Number of particles
particles.nb=particles.Lambda*prod(particles.deltar);

% Generation of the initial positions
particles.R = rand(scheduler.dimr,particles.nb); 
for d=1:scheduler.dimr,
	particles.R(d,:) = particles.deltar(d)*particles.R(d,:)+particles.rmin(d);
end

% Generation of the initial velocities
particles.V=randn(scheduler.dimv,particles.nb);
for d=1:scheduler.dimv,
	particles.V(d,:)=particles.vth(d)*particles.V(d,:);
end

if 0
[N,X]=hist(particles.V(d,:),30);
norm=d01gaf(X(:),N(:));
N=N/norm;
subplot(212)
plot(X,N,X,1/sqrt(2*pi)/particles.vth(d)*exp(-X.^2/2/(particles.vth(d)).^2));
keyboard
end

