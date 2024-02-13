function particles=beam_bnd(particles,scheduler)
% function particles=beam_bnd(particles,dt)

%
% $Id: beam_bnd.m,v 1.5 2011/03/26 15:36:05 patrick Exp $
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

DIMV=[1 2 3];
DIMR=[1:scheduler.dimr];
if particles.nb~=0,

	for d=1:scheduler.dimr,

		ii=find(particles.R(d,:)<=particles.rmin(d));
		nb=length(ii);
		particles.R(:,ii)=[];
		particles.V(:,ii)=[];
		particles.nb=particles.nb-nb;

		ii=find(particles.R(d,:)>=particles.rmax(d));
		nb=length(ii);
		particles.R(:,ii)=[];
		particles.V(:,ii)=[];
		particles.nb=particles.nb-nb;
	end

end

v=zeros(scheduler.dimv,particles.ninject);
r=zeros(scheduler.dimr,particles.ninject);

if particles.u(particles.beamdir) > 0.0,
	v(particles.beamdir,:) = randflux( ...
		particles.u(particles.beamdir),particles.vth(particles.beamdir), ...
		1,particles.ninject);
	r(particles.beamdir,:)=particles.rmin(particles.beamdir)+ ...
		rand(1,particles.ninject).*v(particles.beamdir,:)*scheduler.dt;
else
	v(particles.beamdir,:) = -randflux( ...
		-particles.u(particles.beamdir),particles.vth(particles.beamdir), ...
		1,particles.ninject);
	r(particles.beamdir,:)=particles.rmax(particles.beamdir)+ ...
		rand(1,particles.ninject).*v(particles.beamdir,:)*scheduler.dt;
end

ii=find(DIMV~=particles.beamdir);
for od=DIMV(find(DIMV~=particles.beamdir)),
	v(od,:)=particles.vth(od)*randn(1,particles.ninject)+particles.u(od);
end

for od=DIMR(find(DIMR~=particles.beamdir)),
	r(od,:)=particles.width(od)*randn(1,particles.ninject)+particles.rc(od);
end
particles.V=[particles.V v];
particles.R=[particles.R r];

particles.nb=particles.nb+particles.ninject;

for d=1:scheduler.dimr,
	ii=find(particles.R(d,:)<=particles.rmin(d));
	nb=length(ii);
	particles.R(:,ii)=[];
	particles.V(:,ii)=[];
	particles.nb=particles.nb-nb;
					  
	ii=find(particles.R(d,:)>=particles.rmax(d));
	nb=length(ii);
	particles.R(:,ii)=[];               
	particles.V(:,ii)=[];
	particles.nb=particles.nb-nb;             
end

