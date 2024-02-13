function f=plot2dfield(filename,fieldname,subdim,nrnc,seeds)
% function f=plot2dfield(filename,fieldname,subdim,nrnc,seeds)

%
% $Id: plotphic2d.m,v 1.3 2011/03/26 15:36:05 patrick Exp $
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

f=readhdf(filename,'potential',subdim,seeds);

g=1e-2;
Te=10;

K=g/Te;
if length(f.dims)==3,
	L=max(f.dims{3})-min(f.dims{3});
	[Y,X]=ndgrid(f.dims{3},f.dims{2});
else
	L=max(f.dims{2})-min(f.dims{2});
	[Y,X]=ndgrid(f.dims{2},f.dims{1});
end
n0=K*L/(1-exp(-K*L));
f0=n0.*exp(-K*Y);
y0=1/K*log(n0);
phi=-K*Te*(Y-y0);
if length(f.dims)==3,
	for j=1:size(f.var,3)
		f.var(:,:,j)=f.var(:,:,j)-phi;
	end
else
	f.var = f.var-phi;
end


clim=[min(f.var(:)), max(f.var(:))];

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

for j=1:size(f.var,3),
	subplot(nr,nc,j+jstart), 
	if smoothflag~=0
		imagesc(f.dims{3}, f.dims{2}, ...
			smooth2(f.dims{3},f.dims{2},squeeze(f.var(:,:,j)),...
	  	smoothflag*size(squeeze(f.var(:,:,j)))),clim), 
	else
		if length(f.dims)==3,
			imagesc(f.dims{2}, f.dims{3}, squeeze(f.var(:,:,j)),clim), 
		else
			imagesc(f.dims{1}, f.dims{2}, f.var,clim),
		end
%		pcolor(x,y,squeeze(fxy(:,:,j))), shading interp, caxis(clim)
%		contourf(x,y,squeeze(fxy(:,:,j))),
	end

	if length(f.dims)==3,
		title(sprintf('%s(%s,%s,%s=%.1f)', ...
			f.varname, f.dimsname{2}, f.dimsname{3}, f.dimsname{1}, f.dims{1}(j)))
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

