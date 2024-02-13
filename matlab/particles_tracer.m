function tracer=particles_tracer(tracer,particles,scheduler,pos)
% function tracer=particles_tracer(particles,tracer,pos)

%
% $Id: particles_tracer.m,v 1.5 2011/03/26 15:36:05 patrick Exp $
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

tracer.X(scheduler.citer+1,:)=particles.R(1,tracer.index);
tracer.Y(scheduler.citer+1,:)=particles.R(2,tracer.index);
if scheduler.dimr==3
	tracer.Z(scheduler.citer+1,:)=particles.R(3,tracer.index);
end

subplot(pos)
%cla
if scheduler.citer>1
	ii=1:scheduler.citer+1;
	switch scheduler.dimr,
		case 2,
			plot(tracer.X(ii,:),tracer.Y(ii,:),'o');
		  set(gca,'xlim',[scheduler.rmin(1) scheduler.rmax(1)]);
		  set(gca,'ylim',[scheduler.rmin(2) scheduler.rmax(2)]);

		case 3,
			plot3(tracer.X(ii,:),tracer.Y(ii,:),tracer.Z(ii,:),'o');
			if 1
				subplot(234), plot(tracer.X(ii,:),tracer.Y(ii,:),'o');
				subplot(235), plot(tracer.X(ii,:),tracer.Z(ii,:),'o');
				subplot(236), plot(tracer.Y(ii,:),tracer.Z(ii,:),'o');
			end
		  set(gca,'xlim',[scheduler.rmin(1) scheduler.rmax(1)]);
		  set(gca,'ylim',[scheduler.rmin(2) scheduler.rmax(2)]);
		  set(gca,'zlim',[scheduler.rmin(3) scheduler.rmax(3)]);
	end
	%set(gca,'box','on');
	drawnow
	pause(1)
end

