function [rho,particles]=particles_to_rho(particles,scheduler)
% function [rho,particles]=particles_to_rho(particles,scheduler)

%
% $Id: particles_to_rho.m,v 1.6 2011/03/26 15:36:05 patrick Exp $
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

switch scheduler.dimr
	case 2, 
		rho=zeros(scheduler.nx,scheduler.ny);
	case 3, 
		rho=zeros(scheduler.nx,scheduler.ny,scheduler.nz);
end

lambda=scheduler.Lambda;

for s=1:length(particles),
	fprintf('\tspecies %d: \t# particles %d\n',s,particles{s}.nb);
	if strcmp(particles{s}.type,'beam')
		particles{s}.t(scheduler.citer+1)=scheduler.ctime;
		particles{s}.n(scheduler.citer+1)=particles{s}.nb;
	end
	Z_n0=particles{s}.charge/(scheduler.Lambda/scheduler.unitvolume);
	rho=rho+Z_n0*species_to_rho(particles{s},scheduler);
end

function rho=species_to_rho(species,scheduler)

nx=scheduler.nx;
ny=scheduler.ny;

switch scheduler.dimr
	case 2, 
		rho=zeros(nx,ny);
	case 3, 
		nz=scheduler.nz;
		rho=zeros(nx,ny,nz);
end

hx=scheduler.deltar(1)/(nx-1);
hy=scheduler.deltar(2)/(ny-1);
if scheduler.dimr==3,
	hz=scheduler.deltar(3)/(nz-1);
end

X=(species.R(1,:)-scheduler.rmin(1))/hx;
Y=(species.R(2,:)-scheduler.rmin(2))/hy;
I=fix(X);
J=fix(Y);
DX=X-I;
DY=Y-J;
if scheduler.dimr==3,
	Z=(species.R(3,:)-scheduler.rmin(3))/hz;
	K=fix(Z);
	DZ=Z-K;
end

switch scheduler.dimr
	case 2,
		for p=1:species.nb
			i=I(p)+1; 
			j=J(p)+1;
			dx=DX(p);
			dy=DY(p);
			if i<1,
			  [i,I(p),species.R(1,p),scheduler.rmin(1)]
			end
			if j<1,
			  [j,J(p),species.R(2,p),scheduler.rmin(2)]
			end
			rho([i,i+1],[j,j+1]) = rho([i,i+1],[j,j+1])+[1-dx,dx]'*[1-dy,dy];
		end
		rho([1 end],:)=2*rho([1 end],:);
		rho(:,[1 end])=2*rho(:,[1 end]);

	case 3,
		for p=1:species.nb
			i=I(p)+1; 
			j=J(p)+1;
			k=K(p)+1;
			dx=DX(p);
			dy=DY(p);
			dz=DZ(p);
			rho([i,i+1],[j,j+1],k) = rho([i,i+1],[j,j+1],k) + ...
				[1-dx,dx]'*[1-dy,dy]*(1-dz);
			rho([i,i+1],[j,j+1],k+1) = rho([i,i+1],[j,j+1],k+1) + ...
				[1-dx,dx]'*[1-dy,dy]*dz;
		end
		rho([1 end],:,:)=2*rho([1 end],:,:);
		rho(:,[1 end],:)=2*rho(:,[1 end],:);
		rho(:,:,[1 end])=2*rho(:,:,[1 end]);

end
