function [phi,solver]=solve_rho_to_phi(rho,solver)
% function phi=solve_rho_to_phi(rho,solver)

%
% $Id: solve_rho_to_phi.m,v 1.4 2011/03/26 15:36:05 patrick Exp $
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

if solver.dim==2,
	[phi,solver]=solver2d(rho,solver);
	phi=phi([2:end-1],[2:end-1]);
else
	[phi,solver]=solver3d(rho,solver);
	phi=phi([2:end-1],[2:end-1],[2:end-1]);
end

