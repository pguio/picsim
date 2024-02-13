function kinetic
% function kinetic

%
% $Id: kinetic.m,v 1.2 2011/03/26 15:36:04 patrick Exp $
%
% Copyright (c) 2001-2011 Patrick Guio <patrick.guio@gmail.com>
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

g=-1e-2;

filename='gravity_check.hdf';
if exist(filename,'file')==2,
	filename=which(filename);
else
	error(sprintf('Cannot find file %s', filename));
end

[tt,t,x,vx,fxvx]=read2dfield(filename,'plane2',1);
h=max(x);
dx2=0.5*(x(2)-x(1));
dvx2=0.5*(vx(2)-vx(1));
imagesc(x+dx2,vx+dvx2,fxvx), axis xy
xlabel('y')
ylabel('v_y')
title(sprintf('%s(x,y,\\tau=%d)', tt, t));
colorbar
hold on
v0=sqrt(abs(2*g*h));
h=plot(x,sqrt(2*g*x+v0^2),x,-sqrt(2*g*x+v0^2));
set(h,'LineWidth',2);
hold off

orient landscape; set(gcf,'PaperOrientation','portrait');

if exist('exportfig'),
	exportfig(gcf,'kinetic','format','eps');
else
	print -deps kinetic;
end

