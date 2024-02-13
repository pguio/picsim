function M=makemovie2d2field(filename,field1,field2,subdim,seeds)
% function M=makemovie2d2dfield(filename,field,subdim,seeds)
% To save as mpeg file run mpgwrite(M,M(1).colormap,filename);

%
% $Id: makemovie2d2field.m,v 1.7 2011/03/26 15:36:05 patrick Exp $
%
% Copyright (c) 2000, 2001 Patrick Guio <patrick.guio@gmail.com>
%
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

if ~exist('seeds','var')
  seeds=[];
end

mnf=[];
mxf=[];
for j=1:n
	f=readhdf(filename,field1,{subdim{1}(j), subdim{2:end},seeds);
	mnf=[mnf min(f.var(:))];
	mxf=[mxf max(f.var(:))];
end
clim1=[min(mnf) max(mxf)];

mnf=[];
mxf=[];
for j=1:n
  f=readhdf(filename,field2,{subdim{1}(j), subdim{2:end},seeds);
  mnf=[mnf min(f.var(:))];
  mxf=[mxf max(f.var(:))];
end
clim2=[min(mnf) max(mxf)];

for j=1:length(subdim{1}),

	clf

	f=readhdf(filename,field1,{subdim{1}(j), subdim{2:end}},seeds);
	fxy=f.var;
	x=f.dims{2};
	y=f.dims{3};

	subplot(211), 
	imagesc(x,y,fxy(:,:),clim1); 
	axis xy
	title(sprintf('%s(x,y,\\tau=%d)', f.varname, f.dims{1}))
	colorbar;

	f=readhdf(filename,field2,{subdim{1}(j), subdim{2:end}},seeds);
	fxy=f.var;
	x=f.dims{2};
	y=f.dims{3};

	subplot(212), 
	imagesc(x,y,fxy(:,:),clim2); 
	axis xy
	title(sprintf('%s(x,y,\\tau=%d)', f.varname, f.dims{1}))
	colorbar;

	drawnow
		
	M(j)=getframe(gcf);

end
movie(M,1,2);

