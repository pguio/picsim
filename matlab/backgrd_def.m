function particles=back_def(scheduler)
% function particles=back_def(scheduler)

%
% $Id: backgrd_def.m,v 1.6 2011/03/26 15:36:05 patrick Exp $
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

particles=struct('type',[],'nb',[],'rmax',[],'rmin',[],'vth',[],...
	'R',[],'V',[]);

% Particle type
particles.type='backgrd';

% Particle mass
particles.mass=1;
% Particle charge
particles.charge=1;

% parameter Lambda = number of particles / Debye sphere
switch scheduler.dimr,
	case 2,
		particles.Lambda=30;
	case 3,
		particles.Lambda=10;
end

% Particle temperature
particles.temp=[1 1 1];

% boundary condition type
% -periodic
% -reflective
% -free_space
particles.bc=free_space;
