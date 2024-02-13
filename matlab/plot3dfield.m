function f=plot3dfield(filename,fieldname,subdim,nrnc,plottype,seeds,clim)
% function f=plot3dfield(filename,fieldname,subdim,nrnc,plottype,seeds,clim)

%
% $Id: plot3dfield.m,v 1.10 2011/03/26 15:36:05 patrick Exp $
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

if ~exist('subdim','var'),
	subdim=[];
end

if ~exist('seeds','var')
  seeds=[];
end

f=readhdf(filename,fieldname,subdim,seeds);

if ~exist('clim','var')
	clim=[min(f.var(:)), max(f.var(:))];
end

if ~exist('nrnc','var'),
	nrnc=[1 1];
end

j=1;
nr=nrnc(1);
nc=nrnc(2);
if length(nrnc)>2, 
  jstart=nrnc(3)-1;
else
  jstart=0;
end;

if ~exist('smoothflag'), 
	smoothflag=0; 
end;

if exist('plottype') & ~ischar(plottype)
	icut=fix(plottype{2});
	plottype=plottype{1};
end;

fxyz=f.var;
times=f.dims{1};
x=f.dims{2};
y=f.dims{3};
z=f.dims{4};

for i=1:size(f.var,1),

	subplot(nr,nc,i+jstart), 

	if smoothflag,
		f.var = smooth3(fxyz,'gaussian',5);
	end

	[Y,X,Z]=meshgrid(y,x,z);
	ym=max(y);
	xm=max(x);
	zm=max(z);

	if ~exist('plottype') | isempty('plottype'),
		slice(Y,X,Z,squeeze(fxyz(j,:,:,:)),[1/2]*ym,[1/2]*xm,[1/2]*zm);
		xlabel('y'); 
		ylabel('x'); 
		zlabel('z'); 
		title(sprintf('%s(x,y,z,\\tau=%d)', f.varname, times(j)))
		axis ij;
		shading flat;
	elseif ~isempty('plottype') & strcmp(plottype,'cutxy'),
		if ~exist('icut'),
			icut=fix(size(fxyz,4)/2);
		end
		imagesc(y,x,squeeze(fxyz(j,:,:,icut))); 
		axis xy
		xlabel('y'); 
		ylabel('x');
		title(sprintf('%s(x,y,z=%d,\\tau=%d)', f.varname, icut, times(j)))
	elseif ~isempty('plottype') & strcmp(plottype,'cutxz'),
		if ~exist('icut'),
			icut=fix(size(fxyz,3)/2);
		end
		imagesc(x,z,squeeze(fxyz(j,:,icut,:))'); 
		axis xy
		xlabel('x'); 
		ylabel('z');
		title(sprintf('%s(x,y=%d,z,\\tau=%d)', f.varname, icut, times(j)))
	elseif ~isempty('plottype') & strcmp(plottype,'cutyz'),
		if ~exist('icut'),
			icut=fix(size(fxyz,2)/2);
		end
		imagesc(y,z,squeeze(fxyz(j,icut,:,:))'); 
		axis xy
		xlabel('y'); 
		ylabel('z');
		title(sprintf('%s(x=%d,y,z,\\tau=%d)', f.varname, icut, times(j)))
	elseif ~isempty('plottype') & strcmp(plottype,'slice1'),
		slice(Y,X,Z,squeeze(fxyz(j,:,:,:)),[1/2]*ym,[1/2]*xm,[1/2]*zm);
		xlabel('y'); 
		ylabel('x'); 
		zlabel('z'); 
		title(sprintf('%s(x,y,z,\\tau=%d)', f.varname, times(j)))
		axis ij;
		shading flat;
	elseif ~isempty('plottype') & strcmp(plottype,'slice2')
		slice(Y,X,Z,squeeze(fxyz(j,:,:,:)),[1/3 2/3]*ym,[1/3 2/3]*xm,[1/3 2/3]*zm);
		xlabel('y'); 
		ylabel('x'); 
		zlabel('z'); 
		title(sprintf('%s(x,y,z,\\tau=%d)', f.varname, times(j)))
		axis ij;
		shading flat;
	elseif ~isempty('plottype') & strcmp(plottype,'pcolorxy')
		if ~exist('icut'),
			icut=fix(size(fxyz,4)/2);
		end
		pcolor(y,x,squeeze(fxyz(j,:,:,icut))), 
		shading interp, 
		caxis(clim)
		xlabel('y'); 
		ylabel('x');
		title(sprintf('%s(x,y,z=%d,\\tau=%d)', f.varname, icut, times(j)))
	elseif ~isempty('plottype') & strcmp(plottype,'pcolorxz')
		if ~exist('icut'),
			icut=fix(size(fxyz,3)/2);
		end
		pcolor(x,z,squeeze(fxyz(j,:,icut,:))'), 
		shading interp, 
		caxis(clim)
		xlabel('x'); 
		ylabel('z');
		title(sprintf('%s(x,y=%d,z,\\tau=%d)', f.varname, icut, times(j)))
	elseif ~isempty('plottype') & strcmp(plottype,'pcoloryz')
		if ~exist('icut'),
			icut=fix(size(fxyz,2)/2);
		end
		pcolor(y,z,squeeze(fxyz(j,icut,:,:))'), 
		shading interp, 
		caxis(clim)
		xlabel('y'); 
		ylabel('z');
		title(sprintf('%s(x=%d,y,z,\\tau=%d)', f.varname, icut, times(j)))
	elseif ~isempty('plottype') & strcmp(plottype,'contourxy')
		if ~exist('icut'),
			icut=fix(size(fxyz,4)/2);
		end
		contourf(y,x,squeeze(fxyz(j,:,:,icut))),
		xlabel('y'); 
		ylabel('x');
		title(sprintf('%s(x,y,z=%d,\\tau=%d)', f.varname, icut, times(j)))
	elseif ~isempty('plottype') & strcmp(plottype,'contourxz')
		if ~exist('icut'),
			icut=fix(size(fxyz,3)/2);
		end
		contourf(x,z,squeeze(fxyz(j,:,icut,:))'),
		xlabel('x'); 
		ylabel('z');
		title(sprintf('%s(x,y=%d,z,\\tau=%d)', f.varname, icut, times(j)))
	elseif ~isempty('plottype') & strcmp(plottype,'contouryz')
		if ~exist('icut'),
			icut=fix(size(fxyz,2)/2);
		end
		contourf(y,z,squeeze(fxyz(j,icut,:,:))'),
		xlabel('y'); 
		ylabel('z');
		title(sprintf('%s(x=%d,y,z,\\tau=%d)', f.varname, icut, times(j)))
	elseif ~isempty('plottype') & strcmp(plottype,'volume')
		[Y,X,Z,Fxyz] = subvolume(y,x,z,squeeze(fxyz(j,:,:,:)),...
			 [NaN NaN x(1) x(fix(end/2)+1) z(1) z(fix(end/2)+1)]);
		m=.05;
		m=.0;
		m=1.5;
		clim=[.1 .6];
		p=patch(isosurface(Y,X,Z,Fxyz,m), 'FaceColor','red','EdgeColor','none');
		p2=patch(isocaps(Y,X,Z,Fxyz,m), 'FaceColor','interp','EdgeColor','none');
		view(3); 
		axis ij
		axis vis3d tight;  
		daspect([1 1 1])
		colormap(jet)
		camlight
		lighting gouraud

	  set(gca,'xlim',[min(Y(:)) max(Y(:))]);
	  set(gca,'ylim',[min(X(:)) max(X(:))]);
	  set(gca,'zlim',[min(Z(:)) max(Z(:))]);

		xlabel('y')
		ylabel('x')
		zlabel('z')
		title(sprintf('%s(x,y,z,\\tau=%d)', f.varname, times(j)))
	end


	if 0
		if floor((j+jstart-1)/nc)==nr-1, 
			xlabel('x');
		else
			set(gca,'XTickLabel',[]);
		end
		if mod(j+jstart-1,nc)==0,
			ylabel('y');
		else
			set(gca,'YTickLabel',[]);
		end
		axis xy
		h=colorbar;
		if mod(j+jstart-1,nc)~=nc-1 & length(tidx)>1,
			set(h,'YTickLabel',[]);
		end
		drawnow
	end
	j=j+1;
	drawnow
end



