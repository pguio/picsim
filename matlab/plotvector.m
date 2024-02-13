function f=plotvector(filename,fieldname,subdim,nrnc,seeds)
% function f=plotvector(filename,fieldname,subdim,nrnc,seeds)

%
% $Id: plotvector.m,v 1.7 2011/03/26 15:36:05 patrick Exp $
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

subplot(nr,nc,j+jstart), 
plot(f.dims{1}, f.var)
title(sprintf('%s(%s)', f.varname, f.dimsname{1}))
if floor((j+jstart-1)/nc)==nr-1, 
	xlabel(f.dimsname{1});
else
	set(gca,'XTickLabel',[]);
end
if mod(j+jstart-1,nc)~=0,
	set(gca,'YTickLabel',[]);
end
axis xy
drawnow
