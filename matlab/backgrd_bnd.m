function particles=backgrd_bnd(particles,scheduler)
% function particles=backgrd_bnd(particles,dt)

%
% $Id: backgrd_bnd.m,v 1.6 2011/03/26 15:36:05 patrick Exp $
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

backgrd_const;

switch particles.bc
	case periodic,
		particles=periodic_bnd(particles,scheduler);
	case reflective
		particles=reflective_bnd(particles,scheduler);
	case free_space,
		particles=free_space_bnd(particles,scheduler);
end
	

function  particles=periodic_bnd(particles,scheduler)

for d=1:scheduler.dimr, % for all spatial dimensions

  % find the particles outside the simulation area
  % lower boundary
	ii = find(particles.R(d,:) < particles.rmin(d));
	particles.R(d,ii) = particles.R(d,ii) + particles.deltar(d);

  % upper boundary
	ii = find(particles.R(d,:) >= particles.rmax(d));
  particles.R(d,ii) = particles.R(d,ii) - particles.deltar(d);
end

function  particles=reflective_bnd(particles,scheduler)

for d=1:scheduler.dimr, % for all spatial dimensions

% find the particles outside the simulation minimum
	ii=find(particles.R(d,:) <= particles.rmin(d));
	particles.R(d,ii)=2.0*particles.rmin(d)-particles.R(d,ii);
	particles.V(d,ii)=-particles.V(d,ii);

% max
	ii=find(particles.R(d,:) >= particles.rmax(d));
	particles.R(d,ii)=2.0*particles.rmax(d)-particles.R(d,ii);
	particles.V(d,ii)=-particles.V(d,ii);

end

function particles=free_space_bnd(particles,scheduler)

DIMR=[1:scheduler.dimr];
DIMV=[1:scheduler.dimv];

for d=1:scheduler.dimr, % for all spatial dimensions

% find the particles outside the simulation minimum
	ii=find(particles.R(d,:) <= particles.rmin(d));
	nb=length(ii);
	particles.R(:,ii)=[];
	particles.V(:,ii)=[];

	fprintf(1,'Dim %d minimum %d lost\n', d, nb);

	r=zeros(scheduler.dimr,nb);
	v=zeros(scheduler.dimv,nb);

% distribution along the spatial dimension
	v(d,:) = randflux(0,particles.vth(d),1,nb);
	r(d,:) = particles.rmin(d)+rand(1,nb).*v(d,:)*scheduler.dt;

% distributions in the orthogonal dimensions
	for od=DIMV(find(DIMV~=d)), 
		v(od,:)=particles.vth(od)*randn(1,nb); 
	end;
	for od=DIMR(find(DIMR~=d)), 
		r(od,:)=particles.deltar(od)*rand(1,nb)+particles.rmin(od);
	end
	particles.R=[particles.R r];
	particles.V=[particles.V v];

% max
	ii=find(particles.R(d,:) >= particles.rmax(d));
	nb=length(ii);
	particles.R(:,ii)=[];
	particles.V(:,ii)=[];

	fprintf(1,'Dim %d maximum %d lost\n', d, nb);

	r=zeros(scheduler.dimr,nb);
	v=zeros(scheduler.dimr,nb);

% distribution along the spatial dimension
	v(d,:) = -randflux(0,particles.vth(d),1,nb);
	r(d,:) = particles.rmax(d)+rand(1,nb).*v(d,:)*scheduler.dt;

% distributions in the orthogonal dimensions
	for od=DIMV(find(DIMV~=d)), 
		v(od,:)=particles.vth(od)*randn(1,nb); 
	end;
	for od=DIMR(find(DIMR~=d)), 
		r(od,:)=particles.deltar(od)*rand(1,nb)+particles.rmin(od);
	end
	particles.R=[particles.R r];
	particles.V=[particles.V v];

end
