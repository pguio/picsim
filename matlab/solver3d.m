function [phi,solver]=solver3d(rho,solver)
% function phi=solve_rho_to_phi(rho,solver)

%
% $Id: solver3d.m,v 1.5 2011/03/26 15:36:05 patrick Exp $
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

ngrid=solver.ngrid;
kcycle=solver.kcycle;
nx=solver.nx;
ny=solver.ny;
nz=solver.nz;

phi=zeros(nx+2,ny+2,nz+2);

% set u, u in PHI{}, RHO} and adjust right hand side
[PHI{ngrid},RHO{ngrid}]=swk3(phi,rho,solver);

if solver.iguess==0, % no initial guess at finest grid level!

% transfer down to all grid levels
	for k=ngrid:-1:2,
		[PHI{k-1},RHO{k-1}]=trsfc3(PHI{k},RHO{k},k-1,solver);
	end
% adjust right hand side at all grid levels in case
% rhs or specified b.c. in phi or gbdy changed
	for k=1:ngrid,
		[PHI{k},RHO{k}]=adjmd3(PHI{k},RHO{k},k,solver);
	end
% execute one full multigrid cycle
	PHI{1}=smooth2d(PHI{1},RHO{1},1,solver);

	for k=2:ngrid
		PHI{k}=prolon3(PHI{k-1},PHI{k},k,solver);
	 	PHI=mg3d(PHI,RHO,k,solver);
	end

else,
% transfer down to all grid levels
	for k=ngrid:-1:2,
		[PHI{k-1},RHO{k-1}]=trsfc3(PHI{k},RHO{k},k-1,solver);
	end
	K=ngrid;
% adjust rhs at finest grid level only
	[PHI{K},RHO{K}]=adjmd3(PHI{K},RHO{K},K,solver);
end;

K=ngrid;
% execute maxcy more multigrid k cycles from finest level
epsk=solver.tolmax*mean(abs(RHO{K}(:)));
d=lop3d(PHI{K},K,solver)-RHO{K};
tk=mean(abs(d(:)))-epsk;
if tk<0.0
	fprintf(1,'tk=%.2e mean |defect|=%.2e std |defect|=%.2e\n', ...
	tk,mean(abs(d(:))),std(abs(d(:))))
else
	for i=1:solver.maxcy,
		fprintf(1,'Entering iter#%3d of mg2d with k=%d\n',i,K);
		PHI=mg3d(PHI,RHO,K,solver);
		d=lop3d(PHI{K},K,solver)-RHO{K};
		tk=mean(abs(d(:)))-epsk;
		if tk<0.0
			fprintf(1,'tk=%.2e mean |defect|=%.2e std |defect|=%.2e\n', ...
				tk,mean(abs(d(:))),std(abs(d(:))))
			break;
		end
	end
end

% save defect mean and standard deviation
solver.mean_defect=mean(abs(d(:)));
solver.std_defect=std(abs(d(:)));

phi=PHI{K};
