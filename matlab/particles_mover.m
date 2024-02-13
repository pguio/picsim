function particles=particles_mover(particles,scheduler)
% function particles=particles_mover(particles,solver,rho,scheduler)

%
% $Id: particles_mover.m,v 1.5 2011/03/26 15:36:05 patrick Exp $
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

if scheduler.citer==0,
	for s=1:length(particles)
		particles{s}=species_deccelerate(particles{s},scheduler);
	end
end

for s=1:length(particles)
	particles{s}=species_mover(particles{s},scheduler);
end

function species=species_mover(species,scheduler)

mu=species.mass/scheduler.mass;
Z=species.charge;

if scheduler.electrostatic
	dex=field_to_particles(species,scheduler,scheduler.dE{1});
	dey=field_to_particles(species,scheduler,scheduler.dE{2});
	if 	scheduler.dimr==2,
		dez=zeros(1,species.nb);
	else
		dez=field_to_particles(species,scheduler,scheduler.dE{3});
	end
	E=Z/mu*[scheduler.E(1)+dex; scheduler.E(2)+dey; scheduler.E(3)+dez];
% Translation
	species.V = species.V + 0.5 * E * scheduler.dt;
end

if scheduler.gravitation
% Translation
	species.V = species.V + 0.5 * scheduler.G * scheduler.dt;
end

% Rotation
if  scheduler.electromagnetic==1,
	unit = ones(1,species.nb);
	T = Z/mu*[ scheduler.T(1)*unit; scheduler.T(2)*unit; scheduler.T(3)*unit];
	S = Z/mu*[ scheduler.S(1)*unit; scheduler.S(2)*unit; scheduler.S(3)*unit];

	v = species.V + cross(species.V,T);
	species.V = species.V + cross(v,S);

elseif scheduler.electromagnetic==2,

	if isa(scheduler.extMagneticField, 'function_handle'),
	  [extMagneticField, msg] = fcnchk(scheduler.extMagneticField);
	else
	  error('extMagneticField not a Function Handle');
	end
  if  scheduler.dimr==2,
	  R = {species.R(1,:),species.R(2,:)};
	  B = extMagneticField(R);
	else
	  R = {species.R(1,:),species.R(2,:),species.R(3,:)};
	  B = extMagneticFIeld(R);
	end

  normB = ones(3,1)*sqrt(sum(B.^2));
	T = -Z/mu*tan(normB*scheduler.dt/2.0).*B./normB;
	S = Z/mu*2.0*T./(1+ones(3,1)*sum(T.^2));

	v = species.V + cross(species.V,T);
	species.V = species.V + cross(v,S);

end
	
% Translation
if scheduler.electrostatic,
	species.V = species.V + 0.5 * E * scheduler.dt;
end
if scheduler.gravitation
% Translation
	species.V = species.V + 0.5 * scheduler.G * scheduler.dt;
end

R=species.R;
% Velocity to Position
species.R = species.R + species.V(1:scheduler.dimr,:) * scheduler.dt;


function species=species_deccelerate(species,scheduler,Ex,Ey)

mu=species.mass/scheduler.mass;
Z=species.charge;

if scheduler.electrostatic
	dex=field_to_particles(species,scheduler,scheduler.dE{1});
	dey=field_to_particles(species,scheduler,scheduler.dE{2});
	if 	scheduler.dimr==2,
		dez=zeros(1,species.nb);
	else
		dez=field_to_particles(species,scheduler,scheduler.dE{3});
	end
	E=Z/mu*[scheduler.E(1)+dex; scheduler.E(2)+dey; scheduler.E(3)+dez];
% Translation
	species.V = species.V - 0.5 * E * scheduler.dt;
end

if scheduler.gravitation
% Translation
	species.V = species.V - 0.5 *  scheduler.G * scheduler.dt;
end

% Rotation
if  scheduler.electromagnetic == 1,
	unit = ones(1,species.nb);
	T = Z/mu*[ scheduler.T(1)*unit; scheduler.T(2)*unit; scheduler.T(3)*unit];

	species.V = species.V - cross(species.V,T);

elseif scheduler.electromagnetic == 2,

	if isa(scheduler.extMagneticField, 'function_handle'),
	  [extMagneticField, msg] = fcnchk(scheduler.extMagneticField);
	else
	  error('extMagneticField not a Function Handle');
	end
  if  scheduler.dimr==2,
	  R = {species.R(1,:),species.R(2,:)};
	  B = extMagneticField(R);
	else
	  R = {species.R(1,:),species.R(2,:),species.R(3,:)};
	  B = extMagneticFIeld(R);
	end

  normB = ones(3,1)*sqrt(sum(B.^2));
	T = -Z/mu*tan(normB*scheduler.dt/2.0).*B./normB;

	species.V = species.V + cross(species.V,T);

end

