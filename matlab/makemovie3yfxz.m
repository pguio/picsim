function M=makemovie3yfxz(f1,f2,f3,field,subdim,seeds)
% function M=makemovie3yfxz(f1,f2,f3,field,subdim,seeds)
% To save as mpeg file run mpgwrite(M,M(1).colormap,filename);

%
% $Id: makemovie3yfxz.m,v 1.3 2011/03/26 15:36:05 patrick Exp $
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

f1=readhdf(f1,field,subdim,seeds);

f2=readhdf(f2,field,subdim,seeds);

f3=readhdf(f3,field,subdim,seeds);

clim = [min([f1.var(:);f2.var(:);f3.var(:)]) ...
	max([f1.var(:);f2.var(:);f3.var(:)])];

for j=1:length(subdim{3}),
	subplot(131), 
	imagesc(f1.dims{2},f1.dims{4},squeeze(f1.var(:,j,:))',clim); 
	axis square;
	title(sprintf('%s(x, y=%d, z, \\tau=%.1f)', ...
		f1.varname,f1.dims{3}(j),f1.dims{1}));
	xlabel(f1.dimsname{2});
	ylabel(f1.dimsname{4});

	subplot(132), 
	imagesc(f2.dims{2},f2.dims{4},squeeze(f2.var(:,i,:))',clim); 
	axis square;
	title(sprintf('%s(x, y=%d, z, \\tau=%.1f)', ...
		f2.varname,f2.dims{3}(j),f2.dims{1}));
	xlabel(f2.dimsname{2});

	subplot(133), 
	imagesc(f3.dims{2},f3.dims{4},squeeze(f3.var(:,i,:))',clim); 
	axis square;
	title(sprintf('%s(x, y=%d, z, \\tau=%.1f)', ...
		f3.varname,f3.dims{3}(j),f3.dims{1}));
	xlabel(f3.dimsname{2});

	M(i)=getframe(gcf);
	i=i+1;
end

