function particles=particles_bnd(particles,scheduler)
% function particles=particles_bnd(particles,scheduler)

%
% $Id: particles_bnd.m,v 1.4 2011/03/26 15:36:05 patrick Exp $
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

for s=1:length(particles)

	switch lower(particles{s}.type)
		case 'beam'
			particles{s}=beam_bnd(particles{s},scheduler);
		case 'backgrd'
			particles{s}=backgrd_bnd(particles{s},scheduler);
	end

end
