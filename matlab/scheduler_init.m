function scheduler=scheduler_init(scheduler)
% function scheduler=scheduler_init(scheduler)

%
% $Id: scheduler_init.m,v 1.4 2011/03/26 15:36:05 patrick Exp $
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

scheduler.deltar=scheduler.rmax-scheduler.rmin;

switch scheduler.dimr
	case 2
		scheduler.unitvolume=prod((scheduler.rmax-scheduler.rmin)./ ...
  		[scheduler.nx-1 scheduler.ny-1]);
	case 3
		scheduler.unitvolume=prod((scheduler.rmax-scheduler.rmin)./ ...
		  [scheduler.nx-1 scheduler.ny-1 scheduler.nz-1]);
end

scheduler.T=-tan(norm(scheduler.B)*scheduler.dt/2.0)* ...
	scheduler.B/norm(scheduler.B);
scheduler.S=2.0*scheduler.T/(1+sum((scheduler.T).^2));

rand('state',scheduler.ustate);
randn('state',scheduler.nstate);

