function F=field_to_particles(species,scheduler,field)
% function F=field_to_particles(species,scheduler,field)

%
% $Id: field_to_particles.m,v 1.5 2011/03/26 15:36:05 patrick Exp $
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

if scheduler.dimr==2,
	[nx,ny]=size(field);
else
	[nx,ny,nz]=size(field);
end

F=zeros(1,species.nb);

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

if scheduler.dimr==2,
	for p=1:species.nb
		i=I(p)+1; 
		j=J(p)+1;
		dx=DX(p);
		dy=DY(p);
		F(p)=sum(sum(field([i,i+1],[j,j+1]).*([1-dx,dx]'*[1-dy,dy])));
	end
else
	for p=1:species.nb
		i=I(p)+1; 
		j=J(p)+1;
		k=K(p)+1;
		dx=DX(p);
		dy=DY(p);
		dz=DZ(p);
		F(p)=sum(sum(field([i,i+1],[j,j+1],k).*([1-dx,dx]'*[1-dy,dy]*(1-dz))))+...
			sum(sum(field([i,i+1],[j,j+1],k+1).*([1-dx,dx]'*[1-dy,dy]*dz)));
	end
end
