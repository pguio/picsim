function M=makemovieyfxz(filename,field,subdim,seeds)
% function M=makemovieyfxz(filename,field,subdim,seeds)
% To save as mpeg file run mpgwrite(M,M(1).colormap,filename);

%
% $Id: makemovieyfxz.m,v 1.4 2011/03/26 15:36:05 patrick Exp $
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

f=readhdf(filename,field,subdim,seeds);

if 0
	for j=1:length(subdim{3}),
		f.var(:,j,:)=log10(abs(fftshift(fft2(squeeze(f.var(:,j,:))))));
	end
end

clim = [min(f.var(:)) max(f.var(:))];

for j=1:length(subdim{3}),
	subplot(111), 
	imagesc(f.dims{2},f.dims{4},squeeze(f.var(:,j,:))',clim);
	axis square;
	title(sprintf('%s(x, y=%d, z, \\tau=%.1f)', ...
		f.varname,f.dims{3}(j),f.dims{1}));
	xlabel(f.dimsname{2});
	ylabel(f.dimsname{4});
	M(i)=getframe(gcf);
	i=i+1;
end





