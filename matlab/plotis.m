function [frq,is]=plotis(filename,plotcmd,subdim,nrnc,seeds)
% function [frq,is]=plotis(filename,plotcmd,subdim,nrnc,seeds)

%
% $Id: plotis.m,v 1.5 2011/03/26 15:36:05 patrick Exp $
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

if ~exist('subdim','var'),
  subdim=[];
end
if ~exist('seeds','var'),
  seeds=[];
end

f=readhdf(filename,'avSkw2',subdim,seeds);
is=f.var;
k=f.dims{2};
frq=f.dims{3};

n=length(k);
n2=fix(length(k)/2);

if ~exist('plotcmd')
	plotcmd='plot';
end

for i=1:n2,
	subplot(n2,2,2*i-1)
	feval(plotcmd,sign(k(i))*frq,is(:,i))
	title(sprintf('%s (k=%.2e)',f.varname, k(i)));

	subplot(n2,2,2*i)
	feval(plotcmd, sign(k(n-i+1))*frq, is(:,n-i+1))
	title(sprintf('%s (k=%.2e)', f.varname , k(n-i+1)));
end

