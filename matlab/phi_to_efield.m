function E=phi_to_efield(phi,scheduler,solver)
% function E=phi_to_efield(phi,solver)

%
% $Id: phi_to_efield.m,v 1.6 2011/03/26 15:36:05 patrick Exp $
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

nxa=solver.nxa;
nxb=solver.nxb;
nyc=solver.nyc;
nyd=solver.nyd;

nx=solver.nx;
ny=solver.ny;

hx=(solver.xb-solver.xa)/(nx-1);
hy=(solver.yd-solver.yc)/(ny-1);

if solver.dim==3
	nz=solver.nz;
	hz=(solver.zf-solver.ze)/(nz-1);
	nze=solver.nze;
	nzf=solver.nzf;
end

scheduler_const;
s = solver_const;

switch scheduler.FlagGrad
%  Momentum conserving d/dx(Phi(i))=(phi(i+1)-phi(i-1))/(x(i+1)-x(i-1))
	case momentum_conserving, 
		imn=1;
		nd=2;
%  Energy conserving d/dx(Phi(i))=(phi(i+1)-phi(i))/(x(i+1)-x(i))
	case energy_conserving,
		imn=0;
		nd=1;
end

switch solver.dim
	case 2,
% E_x=-d\phi/dx
		I=[2:nx-1];
		J=[1:ny];
		Ex=zeros(nx,ny);
		Ex(I,:)=-(phi(I+1,J)-phi(I-imn,J))/(nd*hx);
% Periodic boundary condition at x=xa and x=xb
		if nxa==bc.periodic & nxb==bc.periodic,
			Ex(1,:)=-(phi(2,J)-phi(nx-imn,J))/(nd*hx);
			Ex(nx,:)=Ex(1,:);
		else
			Ex(1,:)=-(phi(2,J)-phi(1,J))/hx;
			Ex(nx,:)=-(phi(nx,J)-phi(nx-1,J))/hx;
		end
% E_y=-d\phi/dy
		I=[1:nx];
		J=[2:ny-1];
		Ey=zeros(nx,ny);
		Ey(:,J)=-(phi(I,J+1)-phi(I,J-imn))/(nd*hy);
% Periodic boundary condition at y=yc and y=yd
		if nyc==bc.periodic & nyd==bc.periodic,
			Ey(:,1)=-(phi(I,2)-phi(I,ny-imn))/(nd*hy);
			Ey(:,ny)=Ey(:,1);
		else
			Ey(:,1)=-(phi(I,2)-phi(I,1))/hy;
			Ey(:,ny)=-(phi(I,ny)-phi(I,ny-1))/hy;
		end
% Form vector E
	E={Ex, Ey};

	case 3,
% E_x=-d\phi/dx
		I=[2:nx-1];
		J=[1:ny];
		K=[1:nz];
		Ex=zeros(nx,ny,nz);
		Ex(I,:,:)=-(phi(I+1,J,K)-phi(I-imn,J,K))/(nd*hx);
% Periodic boundary condition at x=xa and x=xb
		if nxa==bc.periodic & nxb==bc.periodic,
			Ex(1,:,:)=-(phi(2,J,K)-phi(nx-imn,J,K))/(nd*hx);
			Ex(nx,:,:)=Ex(1,:,:);
		else
			Ex(1,:,:)=-(phi(2,J,K)-phi(1,J,K))/hx;
			Ex(nx,:,:)=-(phi(nx,J,K)-phi(nx-1,J,K))/hx;
		end
% E_y=-d\phi/dy
		I=[1:nx];
		J=[2:ny-1];
		Ey=zeros(nx,ny,nz);
		Ey(:,J,:)=-(phi(I,J+1,K)-phi(I,J-imn,K))/(nd*hy);
% Periodic boundary condition at y=yc and y=yd
		if nyc==bc.periodic & nyd==bc.periodic,
			Ey(:,1,:)=-(phi(I,2,K)-phi(I,ny-imn,K))/(nd*hy);
			Ey(:,ny,:)=Ey(:,1,:);
		else
			Ey(:,1,:)=-(phi(I,2,K)-phi(I,1,K))/hy;
			Ey(:,ny,:)=-(phi(I,ny,K)-phi(I,ny-1,K))/hy;
		end
% E_z=-d\phi/dz
		J=[1:ny];
		K=[2:nz-1];
		Ez=zeros(nx,ny,nz);
		Ez(:,:,K)=-(phi(I,J,K+1)-phi(I,J,K-imn))/(nd*hz);
% Periodic boundary condition at z=ze and z=zf
		if nze==bc.periodic & nzf==bc.periodic,
			Ez(:,:,1)=-(phi(I,J,2)-phi(I,J,nz-imn))/(nd*hz);
			Ez(:,:,nz)=Ez(:,:,1);
		else
			Ez(:,:,1)=-(phi(I,J,2)-phi(I,J,1))/hz;
			Ez(:,:,nz)=-(phi(I,J,nz)-phi(I,J,nz-1))/hz;
		end
% Form vector E
		E={Ex, Ey, Ez};
end

