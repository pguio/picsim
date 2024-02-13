function M=makemovie2d(filename,field,subdim,seeds,clim)
% function M=makemovie2d(filename,field,subdim,seeds,clim)
% To save as mpeg file run mpgwrite(M,M(1).colormap,filename);

%
% $Id: makemovie2d.m,v 1.10 2011/03/26 15:36:05 patrick Exp $
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

if ~exist('seeds','var')
	seeds=[];
end

mnf=[];
mxf=[];
if ~exist('clim','var'),
	for j=1:length(subdim{1}),
		f=readhdf(filename,field,{subdim{1}(j), subdim{2:end}},seeds);
		%	[tt,times,x,y,fxy]=meanhdf2dfield(filename,seeds,field,j);

		mnf=[mnf min(f.var(:))];
		mxf=[mxf max(f.var(:))];
	end
	clim=[min(mnf) max(mxf)];
end

for j=1:length(subdim{1}),

	f=readhdf(filename,field,{subdim{1}(j), subdim{2:end}},seeds);
	fxy=f.var;
	x=f.dims{2};
	y=f.dims{3};
	clf
	imagesc(x,y,squeeze(fxy(1,:,:))',clim); 
	axis xy
	axis equal
	axis tight
	title(sprintf('%s(x,y,\\tau=%.3f)', f.varname, f.dims{1}))
	colorbar;
	drawnow
		
	M(j)=getframe(gcf);

end
%movie(M,1,2);

