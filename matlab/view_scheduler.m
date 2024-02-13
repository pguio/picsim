function view_scheduler(scheduler)
% function view_scheduler(scheduler)

%
% $Id: view_scheduler.m,v 1.4 2011/03/26 15:36:06 patrick Exp $
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

fprintf(1,'\n************************\n');
fprintf(1,'DYNAMIC SCHEDULER SET UP\n');
fprintf(1,'************************\n');
fprintf(1,'%d SPATIAL DIMENSIONS\n', scheduler.dimr);
fprintf(1,'%d VELOCITY DIMENSIONS\n', scheduler.dimv);

fprintf(1,'UNIFORM RANDOM STATE=%d\n', scheduler.ustate);
fprintf(1,'NORMAL  RANDOM STATE=%d\n', scheduler.nstate);

fprintf(1,'TIME INCREMENT= %.1f\n', scheduler.dt);
fprintf(1,'MAX TIME STEPS= %d\n', scheduler.maxiter);

fprintf(1,'NORMALISATION MASS   =%.2f\n', scheduler.mass);
fprintf(1,'NORMALISATION LAMBDA =%.2f\n', scheduler.Lambda);
fprintf(1,'NORMALISATION TEMP   =%.2f\n', scheduler.temp);

fprintf(1,'NX=%3d ==> XMIN=%5.1f, XMAX=%5.1f\n', ...
	scheduler.nx, scheduler.rmin(1), scheduler.rmax(1));
fprintf(1,'NY=%3d ==> YMIN=%5.1f, YMAX=%5.1f\n', ...
	scheduler.ny, scheduler.rmin(2), scheduler.rmax(2));
if scheduler.dimr==3,
	fprintf(1,'NZ=%3d ==> ZMIN=%5.1f, ZMAX=%5.1f\n', ...
		scheduler.nz, scheduler.rmin(3), scheduler.rmax(3));
end

msg_electrostatic(scheduler);
msg_electromagnetic(scheduler);
msg_gravitation(scheduler);

fprintf(1,'****************\n');

function msg_electrostatic(scheduler)

scheduler_const;

if scheduler.electrostatic,
	fprintf(1,'ELECTROSTATIC FORCE   ON\n');
	fprintf(1,'==> GRADIENT = %s\n', FlagGradType{scheduler.FlagGrad});
else
	fprintf(1,'ELECTROSTATIC FORCE   OFF\n');
end

function msg_electromagnetic(scheduler)

if scheduler.electromagnetic,
	fprintf(1,'ELECTROMAGNETIC FORCE ON\n');
else
	fprintf(1,'ELECTROMAGNETIC FORCE OFF\n');
end
function msg_gravitation(scheduler)

if scheduler.gravitation,
	fprintf(1,'GRAVITATION FORCE     ON\n');
else
	fprintf(1,'GRAVITATION FORCE     OFF\n');
end
