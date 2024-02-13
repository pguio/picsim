function particles=beam_def(scheduler)
% function particles=beam_def

%
% $Id: beam_def.m,v 1.4 2011/03/26 15:36:05 patrick Exp $
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

particles=struct('nb',[],'rmax',[],'rmin',[],...
	'u',[],'vth',[],...
  'R',[],'V',[]);

particles.type='beam';

% Particle mass
particles.mass=1;
% Particle charge
particles.charge=1;

switch scheduler.dimr
	case 2
		particles.Lambda=10;
	case 3
		particles.Lambda=2;
end

particles.beamdir=1;

particles.percent=0.988;

particles.u=-[4 4 4];
particles.temp=[0.5 0.5 0.5];

particles.rc=0.5*(scheduler.rmin+scheduler.rmax);
switch scheduler.dimr
	case 2
		particles.seff=[20 20];
		particles.depth=[40 40];
	case 3
		particles.seff=[5 5 5];
		particles.depth=[10 10 10];
end


