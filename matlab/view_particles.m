function view_particles(particles,scheduler)
% function view_particles(particles)

%
% $Id: view_particles.m,v 1.4 2011/03/26 15:36:06 patrick Exp $
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

fprintf(1,'\n******************\n');
fprintf(1,'SPECIES DEFINITION\n');
fprintf(1,'******************\n');
for s=1:length(particles)
	fprintf(1,'\nSPECIES %2d ==> %s\n', s, upper(particles{s}.type));
	fprintf(1,'----------\n');
	switch particles{s}.type
		case 'backgrd',
			view_backgrd(particles{s},scheduler)
		case 'beam',
			view_beam(particles{s},scheduler)
	end
	fprintf(1,'\n');
end
fprintf(1,'******************\n');

function view_backgrd(species,scheduler)

backgrd_const;

fprintf(1,'MASS  =%.2f\n', species.mass);
fprintf(1,'CHARGE=%.2f\n', species.charge);
fprintf(1,'LAMBDA=%d ==> PARTICLE NUMBER=%d\n', species.Lambda, species.nb);

fprintf(1,'ALONG X: TEMP=%5.1f\n', species.temp(1));
fprintf(1,'ALONG Y: TEMP=%5.1f\n', species.temp(2));
fprintf(1,'ALONG Z: TEMP=%5.1f\n', species.temp(3));

fprintf(1,'BOUNDARY TYPE= %s\n',BCTypeName{species.bc+1});


function view_beam(species,scheduler)

AXISNAME={'X','Y','Z'};
DIMR=[1:scheduler.dimr];

fprintf(1,'MASS  =%.2f\n', species.mass);
fprintf(1,'CHARGE=%.2f\n', species.charge);

fprintf(1,'BEAM DIRECTION: %s\n', AXISNAME{species.beamdir});
fprintf(1,'FLUX=%5.1f\n', species.flux(species.beamdir));
switch scheduler.dimr
	case 2,
		fprintf(1,'SECTION=%5.1f\n', species.seff(species.beamdir));
	case 3,
		fprintf(1,'SECTION(%s x %s)=%5.1f x%5.1f\n', ...
			AXISNAME{DIMR(find(DIMR~=species.beamdir))}, ...
			species.seff(DIMR(find(DIMR~=species.beamdir))));
end
fprintf(1,'LAMBDA=%.1f ==>\n', species.Lambda);
fprintf(1,'\tNINIT=%d\n\tNINJECT=%d\n', species.ninit, species.ninject);


fprintf(1,'ALONG X: DRIFT=%5.1f, TEMP=%5.1f\n', species.u(1), species.temp(1));
fprintf(1,'ALONG Y: DRIFT=%5.1f, TEMP=%5.1f\n', species.u(2), species.temp(2));
fprintf(1,'ALONG Z: DRIFT=%5.1f, TEMP=%5.1f\n', species.u(3), species.temp(3));
