function M=makemovieslicer(filename,n)
% function M=makemovieslicer(filename,n)
% To save as mpeg file run mpgwrite(M,M(1).colormap,filename);

%
% $Id: makemovieslicer.m,v 1.3 2011/03/26 15:36:05 patrick Exp $
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

mnf=[];
mxf=[];
for j=1:n
	[tt,dim,subf,parameters]=readslicer(filename,j);
	% [tt,times,x,y,fxy]=mean2dfield(filename,seeds,field,j);
	mnf=[mnf min(subf(:))];
	mxf=[mxf max(subf(:))];
end
clim=[min(mnf) max(mxf)];

xlim = parameters(1:2);
ylim = parameters(3:4);
zlim = parameters(5:6);

vxlim = parameters(7:8);
vylim = parameters(9:10);
vzlim = parameters(11:12);

if diff(vxlim)==0 & diff(vylim)==0,
  Xlim = zlim;
	Ylim = vzlim;
	Xlabel='z';
	Ylabel='v_z';
elseif diff(vxlim)==0 & diff(vzlim)==0,
	Xlim = ylim;
	Ylim = vylim;
	Xlabel='y';
	Ylabel='v_y';
elseif diff(vylim)==0 & diff(vzlim)==0,
	Xlim = xlim;
	Ylim = vxlim;
	Xlabel='x';
	Ylabel='v_x';
end
x = linspace(Xlim(1), Xlim(2), dim(1));
y = linspace(Ylim(1), Ylim(2), dim(2));

for j=1:n
	[tt,dim,subf,parameters]=readslicer(filename,j);
	% [tt,times,x,y,fxy]=mean2dfield(filename,seeds,field,j);
	clf
	imagesc(x,y,subf), axis xy; 
	set(gca,'clim',clim);
	title(sprintf('time=%d', parameters(13)));
	colormap(hsv)
	%colorbar;
	drawnow
	M(j)=getframe(gcf);
end
movie(M,1,2);


