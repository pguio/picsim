function particles_plot(particles,scheduler,pos)
% function particles_plot(particles,scheduler,pos)

%
% $Id: particles_plot.m,v 1.5 2011/03/26 15:36:05 patrick Exp $
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

color=[
[0 0 1]; ...
[0 .5 0]; ...
[1 0 0]; ...
[0 .75 .75]; ...
[.75 0  .75]; ...
[.75 .75 0]; ...
[.25, .25 ,.25]; ...
];

clf
subplot(pos),
%cla
%view(scheduler.dimr)
%set(gca,'box','on');
hold on;
switch scheduler.dimr,
	case 2,
		for s=1:length(particles),
			if particles{s}.nb~=0,
				h=plot(particles{s}.R(1,:),particles{s}.R(2,:),'.');
				set(h,'color', color(mod(s-1,7)+1,:));
			end
		end
		set(gca,'xlim',[scheduler.rmin(1) scheduler.rmax(1)]);
		set(gca,'ylim',[scheduler.rmin(2) scheduler.rmax(2)]);
	
	case 3,
		for s=1:length(particles),
			if particles{s}.nb~=0,
				h=plot3(particles{s}.R(1,:),particles{s}.R(2,:),...
					particles{s}.R(3,:),'.');
				set(h,'color', color(mod(s-1,7)+1,:));
			end
		end
		set(gca,'xlim',[scheduler.rmin(1) scheduler.rmax(1)]);
		set(gca,'ylim',[scheduler.rmin(2) scheduler.rmax(2)]);
		set(gca,'zlim',[scheduler.rmin(3) scheduler.rmax(3)]);
end
title(sprintf('Particles iter=%d t=%.1f',scheduler.citer,scheduler.ctime));
hold off;

drawnow
