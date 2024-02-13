function f=plot2dfield(filename,fieldname,subdim,nrnc,seeds,clim)
% function f=plot2dfield(filename,fieldname,subdim,nrnc,seeds,clim)

%
% $Id: plot2dfield.m,v 1.12 2011/03/26 15:36:05 patrick Exp $
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
if ~exist('seeds','var'),
	seeds=[];
end

if ~ishdf(filename),
	f=readhdf(filename,fieldname,subdim,seeds);
else
	f=filename;
end

if ~exist('clim','var'),
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

if ~exist('smoothflag','var'), 
	smoothflag=0; 
end;

% time is always first dimension
if ndims(f.var)>2,
	nj=size(f.var,1);
else
	nj=1;
end

for j=1:nj,
	subplot(nr,nc,j+jstart), 
	if smoothflag~=0
		imagesc(f.dims{3}, f.dims{2}, ...
			smooth2(f.dims{3},f.dims{2},squeeze(f.var(:,:,j)),...
	  	smoothflag*size(squeeze(fxy(:,:,j)))),clim), 
	else
		if length(f.dims)==4,
			switch ndims(f.var),
				case 2, imagesc(f.dims{1}, f.dims{2}, f.var, clim), 
				case 3, imagesc(f.dims{2}, f.dims{3}, squeeze(f.var(j,:,:))', clim), 
			end
		else
			switch ndims(f.var),
				case 2, imagesc(f.dims{1}, f.dims{2}, f.var', clim),
				case 3, imagesc(f.dims{2}, f.dims{3}, squeeze(f.var(j,:,:))', clim),
			end
		end
%		pcolor(x,y,squeeze(fxy(:,:,j))), shading interp, caxis(clim)
%		contourf(x,y,squeeze(fxy(:,:,j))),
	end

	if length(f.dims)==4 | ndims(f.var)>2,
		title(sprintf('%s(%s,%s,%s=%.1f)', f.varname, ...
			f.dimsname{2}, f.dimsname{3},  f.dimsname{1}, f.dims{1}(j)))
	else
		title(sprintf('%s(%s,%s)', f.varname, f.dimsname{1},f.dimsname{2}))
	end
	if floor((j+jstart-1)/nc)==nr-1, 
		if length(f.dims)==3,
			xlabel(f.dimsname{2});
		else
			xlabel(f.dimsname{1});
		end
	else
		set(gca,'XTickLabel',[]);
	end
	if mod(j+jstart-1,nc)==0,
		if length(f.dims)==3,
			ylabel(f.dimsname{3});
		else
			ylabel(f.dimsname{2});
		end
	else
		set(gca,'YTickLabel',[]);
	end
	axis xy
	h=colorbar;
	if mod(j+jstart-1,nc)~=nc-1 & size(f.var,3)>1,
		set(h,'YTickLabel',[]);
	end
	drawnow
	j=j+1;
end

