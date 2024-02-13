function particles=beam_init(particles,scheduler)
% function particles=beam_init(particles,scheduler)

%
% $Id: beam_init.m,v 1.6 2011/03/26 15:36:05 patrick Exp $
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

DIMR=[1:scheduler.dimr];
DIMV=[1:scheduler.dimv];

particles.vth=sqrt(particles.temp/particles.mass);

particles.u(DIMR(DIMR~=particles.beamdir))=0.0;

particles.flux=particles.Lambda*particles.u;

particles.width=particles.seff/sqrt(2)/erfinv(particles.percent);

particles.T=particles.depth(particles.beamdir)/ ...
	abs(particles.u(particles.beamdir));

od=DIMR(DIMR~=particles.beamdir);
switch scheduler.dimr,
	case 2, seff=particles.seff(od);
	case 3, seff=pi/4*prod(particles.seff(od));
end

particles.ninit=fix(abs(particles.flux(particles.beamdir))*seff*particles.T);
particles.ninject=fix(abs(particles.flux(particles.beamdir))*seff*scheduler.dt);

particles.rmin=scheduler.rmin;
particles.rmax=scheduler.rmax;

if particles.ninit~=0,

	particles.R=zeros(scheduler.dimr,particles.ninit);
	particles.V=zeros(scheduler.dimv,particles.ninit);

	if particles.u(particles.beamdir) > 0.0,
		particles.V(particles.beamdir,:) = randflux( ...
			particles.u(particles.beamdir),particles.vth(particles.beamdir),...
			1,particles.ninit);
		particles.R(particles.beamdir,:)=particles.rmin(particles.beamdir)+ ...
			rand(1,particles.ninit).*particles.V(particles.beamdir,:)*particles.T;
	else
		particles.V(particles.beamdir,:) = -randflux( ...
			-particles.u(particles.beamdir),particles.vth(particles.beamdir),...
			1,particles.ninit);
		particles.R(particles.beamdir,:)=particles.rmax(particles.beamdir)+ ...
			rand(1,particles.ninit).*particles.V(particles.beamdir,:)*particles.T;
	end

	for od=DIMV(find(DIMV~=particles.beamdir)),
		particles.V(od,:)=particles.vth(od)*randn(1,particles.ninit)+...
			particles.u(od);
	end

	for od=DIMR(find(DIMR~=particles.beamdir)),
		particles.R(od,:)=particles.width(od)*randn(1,particles.ninit)+...
			particles.rc(od);
	end;

end

particles.nb=particles.ninit;

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

