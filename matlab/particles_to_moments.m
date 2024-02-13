function [moments,particles]=particles_to_moments(particles,scheduler)
% function [rho,particles]=particles_to_moments(particles,scheduler)

%
% $Id: particles_to_moments.m,v 1.4 2011/03/26 15:36:05 patrick Exp $
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

nx=scheduler.nx;
ny=scheduler.ny;

switch scheduler.dimr
	case 2, 
		pnumber=zeros(nx,ny);
		for m=scheduler.moments;
			moments{m}=zeros(2,nx,ny);
		end
	case 3, 
		nz=scheduler.nz;
		pnumber=zeros(nx,ny,nz);
		for m=scheduler.moments;
			moments{m}=zeros(3,nx,ny,nz);
		end
end

for s=1:length(particles),
	[snumber,smoments]=species_to_moments(particles{s},scheduler);
	pnumber=pnumber+snumber;
	for m=scheduler.moments
		moments{m}=moments{m}+smoments{m};
	end
end

moments=normalise_moments(pnumber,moments,scheduler);

function moments=normalise_moments(pnumber,moments,scheduler)

pnumber(find(pnumber<=1))=Inf;

switch scheduler.dimr
	case 2, 
		for d=1:scheduler.dimr,
			moments{1}(d,:,:)=squeeze(moments{1}(d,:,:))./pnumber(:,:);
			moments{2}(d,:,:)=squeeze(moments{2}(d,:,:))./(pnumber(:,:)-1);
		end
	case 3,
		for d=1:scheduler.dimr,
			moments{1}(d,:,:,:)=squeeze(moments{1}(d,:,:,:))./pnumber(:,:,:);
			moments{2}(d,:,:,:)=squeeze(moments{2}(d,:,:,:))./(pnumber(:,:,:)-1);
		end
end

function [pnumber,moments]=species_to_moments(species,scheduler)

nx=scheduler.nx;
ny=scheduler.ny;

switch scheduler.dimr
	case 2, 
		for m=scheduler.moments,
			moments{m}=zeros(2,nx,ny);
		end
		pnumber=zeros(nx,ny);
	case 3, 
		nz=scheduler.nz-1;
		for m=scheduler.moments,
			moments{m}=zeros(3,nx,ny,nz);
		end
		pnumber=zeros(nx,ny,nz);
end

hx=scheduler.deltar(1)/(nx-1);
hy=scheduler.deltar(2)/(ny-1);
if scheduler.dimr==3,
	hz=scheduler.deltar(3)/(nz-1);
end

X=(species.R(1,:)-scheduler.rmin(1))/hx;
Y=(species.R(2,:)-scheduler.rmin(2))/hy;
I=fix(X+0.5);
J=fix(Y+0.5);
if scheduler.dimr==3,
	Z=(species.R(3,:)-scheduler.rmin(3))/hz;
	K=fix(Z+0.5);
end

switch scheduler.dimr
	case 2,
		for p=1:species.nb
			i=I(p)+1; 
			j=J(p)+1;
			pnumber(i,j)=pnumber(i,j)+1;
			for m=scheduler.moments,
				moments{m}(1,i,j)=moments{m}(1,i,j)+power(species.V(1,p),m);
				moments{m}(2,i,j)=moments{m}(2,i,j)+power(species.V(2,p),m);
			end
		end

	case 3,
		for p=1:species.nb
			i=I(p)+1; 
			j=J(p)+1;
			k=K(p)+1;
			pnumber(i,j,k)=pnumber(i,j,k)+1;
			for m=scheduler.moments,
				moments{m}(1,i,j,k)=moments{m}(1,i,j,k)+power(species.V(1,p),m);
				moments{m}(2,i,j,k)=moments{m}(2,i,j,k)+power(species.V(2,p),m);
				moments{m}(3,i,j,k)=moments{m}(3,i,j,k)+power(species.V(3,p),m);
			end
		end
	
end

