function M=makemovie3d(filename,field,subdim,seeds)
% function M=makemovie3d(filename,field,subdim,seeds)
% To save as mpeg file run mpgwrite(M,M(1).colormap,filename);

%
% $Id: makemovie3d.m,v 1.5 2011/03/26 15:36:05 patrick Exp $
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

if ~exist('seeds','var')
  seeds=[];
end

mnf=[];
mxf=[];
meanf=[];
for j=1:length(subdim{1}),
	f=readhdf(filename,field,{subdim{1}(j), subdim{2:end},seeds);

	mnf=[mnf min(f.var(:))];
	mxf=[mxf max(f.var(:))];
	meanf=[meanf mean(f.var(:))];
	fprintf(1,'Image %d clim=[%.2f, %.2f]\n', j, min(f.var(:)), max(f.var(:)));
end
clim=[min(mnf) max(mxf)];

%m=mean(meanf);
m=.05;
%clim=[.05 .55];
clim=[.05 .65];
m=.25;
clim=[.25 .55];
m=.15;
clim=[.15 .3];

plottype=[];

for j=1:length(subdim{1}),

	f=readhdf(filename,field,{subdim{1}(j), subdim{2:end}},seeds);
	fxy=f.var;
	x=f.dims{2}(:);
	y=f.dims{3}(:);
	z=f.dims{4}(:);
	clf
	[X,Y,Z]=meshgrid(y,x,z);
	ym=max(y);
	xm=max(x);
	zm=max(z));

	if ~exist('plottype') | isempty('plottype'),
		imagesc(fxyz(:,:,fix(end/2))'); 
	elseif ~isempty('plottype') & strcmp(plottype,'slice2')
		slice(X,Y,Z,fxyz,[1/3 2/3]*ym,[1/3 2/3]*xm,[1/3 2/3]*zm);
	elseif ~isempty('plottype') & strcmp(plottype,'slice1')
		slice(X,Y,Z,fxyz,[]*ym,[1/2]*xm,[1/2]*zm); shading interp
	elseif ~isempty('plottype') & strcmp(plottype,'pcolor')
		pcolor(fxyz(:,:,fix(end/2))'); 
		shading interp
	set(gca,'clim',clim);
	elseif ~isempty('plottype') & strcmp(plottype,'volume')
		%fxyz(ind2sub(size(fxyz),find(X<x(fix(end/2)+1) & Z<z(fix(end/2)+1))))=NaN;
		[X,Y,Z,fxyz] = subvolume(y,x,z,fxyz,...
	 		[NaN NaN x(fix(end/2)+1) x(end) z(1) z(fix(end/2)+1)]);
		fxyz(find(fxyz>clim(2)))=NaN;
		p=patch(isosurface(X,Y,Z,fxyz,m),...
			'FaceColor','red','EdgeColor','none');
		p2=patch(isocaps(X,Y,Z,fxyz,m),...
			'FaceColor','interp','EdgeColor','none');
		set(gca,'xlim',[min(X(:)) max(X(:))]);
		set(gca,'ylim',[min(Y(:)) max(Y(:))]);
		set(gca,'zlim',[min(Z(:)) max(Z(:))]);
		set(gca,'clim',clim);
		view(3); axis vis3d %tight;  
		daspect([1 1 1])
		colormap(hsv)
		camlight
		lighting gouraud

		colorbar;
	end
		
	M(j)=getframe(gcf);
	fprintf(1,'Image %d (%dx%dx%d)\n', j, size(M(j).cdata));

end
movie(M,1,2);

% mpgwrite(M,colormap(jet),mpgfilename)
