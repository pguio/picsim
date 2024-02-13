function [particles,scheduler]=picsim3d(initfile,varargin)
% function [particles,scheduler]=picsim3d(initfile,varargin)

%
% $Id: picsim3d.m,v 1.8 2011/03/26 15:36:05 patrick Exp $
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

close all

[scheduler,particles] = feval(['picsim_' initfile '-3d'], varargin{:});
particles_plot(particles,scheduler,111);

view_scheduler(scheduler);
view_particles(particles,scheduler);

solver=solver_init(scheduler);
if ~isempty(solver),
	solver=solver_par_parsing(solver,varargin{:});
	view_solver(solver);
end

if scheduler.tracer,
  tracer=tracer_init(scheduler);
end

if 0
[rho,particles]=particles_to_rho(particles,scheduler);
field_plot(rho,'\rho',scheduler,111);
[moments,particles]=particles_to_moments(particles,scheduler);
keyboard
end

for i=1:scheduler.maxiter,

	s=sprintf('iter %d/%d time %.1f\n', ...
		scheduler.citer,scheduler.maxiter,scheduler.ctime);
	fprintf(1,'%s%s', s, char(repmat(8,size(s))));

	[rho,particles]=particles_to_rho(particles,scheduler);
	field_plot(rho,'\rho',scheduler,212);
	if scheduler.electrostatic,
		 scheduler=rho_to_efield(rho,scheduler,solver);
	end

	particles=particles_mover(particles,scheduler);

	particles=particles_bnd(particles,scheduler);

	particles_plot(particles,scheduler,111);
	if scheduler.tracer,
	  tracer=particles_tracer(tracer,particles{1},scheduler,212);
	end

	scheduler=scheduler_update(scheduler);

end

if scheduler.maxiter,
	fprintf(1,'iter %d/%d time %.1f\n', ...
		scheduler.citer,scheduler.maxiter,scheduler.ctime);
end

for s=1:length(particles)
	if strcmp(particles{s}.type,'beam')
		beam_stat(particles{s},scheduler);
	end
end
