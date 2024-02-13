function solver=solver_init(scheduler)
% function solver=solver_init(scheduler)

%
% $Id: solver_init.m,v 1.6 2011/03/26 15:36:05 patrick Exp $
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

s = solver_const;

bc = s.bc.mixed;

if scheduler.electrostatic,

	switch scheduler.dimr

		case 2
  		solver=solver2d_def;
  		solver=solver_par_parsing(solver,'nx',scheduler.nx,'ny',scheduler.ny,...
	    	'xa',scheduler.rmin(1),'xb',scheduler.rmax(1),...
	    	'yc',scheduler.rmin(2),'yd',scheduler.rmax(2),...
				'nxa',bc,'nxb',bc,'nyc',bc,'nyd',bc,...
	    	'temin',2,'temax',2,'display',-1);
	  	solver=initpde2(solver);

		case 3
			solver=solver3d_def;
  		solver=solver_par_parsing(solver,...
				'nx',scheduler.nx,'ny',scheduler.ny,'nz',scheduler.nz,...
	    	'xa',scheduler.rmin(1),'xb',scheduler.rmax(1),...
	    	'yc',scheduler.rmin(2),'yd',scheduler.rmax(2),...
	    	'ze',scheduler.rmin(3),'zf',scheduler.rmax(3),...
	    	'temin',2,'temax',2,'display',-1);
			solver=initpde3(solver);
	end

else

  solver=[];

end

