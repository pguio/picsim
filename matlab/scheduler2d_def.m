function scheduler=scheduler2d_def
% function scheduler=scheduler2d_def

%
% $Id: scheduler2d_def.m,v 1.5 2011/03/26 15:36:05 patrick Exp $
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

scheduler_const;

% Dimensionality 
scheduler.dimr=2; % space
scheduler.dimv=3; % velocity

% Uniform random generator initial state
scheduler.ustate=0;
% Normal random generator initial state
scheduler.nstate=0;

% Time increment
scheduler.dt=0.2;
% Max number of iterations
scheduler.maxiter=50;
% Current iteration 
scheduler.citer=0;
% Current time
scheduler.ctime=0;

% Normalisation mass, 
scheduler.mass=1;
% Normalisation plasma parameter, 
scheduler.Lambda=20;
% Normalisation temperature, 
scheduler.temp=1;

% Simulation domain limit
scheduler.rmin=[0 0];
scheduler.rmax=[80 112];

% Number of gridpoints (to calculate rho, phi, E, etc)
scheduler.nx=81;
scheduler.ny=113;

% Moments to calculate
scheduler.moments=[1 2];

% Tracing particles trajectory
scheduler.tracer = false;

% Electrostatic force
scheduler.electrostatic=1;
% Vector of extern electrostatic field
scheduler.E=[0.0 0.0 0.0];
% Calculation of E=grad Phi
scheduler.FlagGrad=momentum_conserving;

% Electromagnetic force
scheduler.electromagnetic=1;
% B Vector 
scheduler.B=[0.0 1.0 0.0];
% external magnetic field function
scheduler.extMagneticField='';

% Gravitation force
scheduler.gravitation=0;
% Vector gravitation
scheduler.G=[0.0 1.0 0.0];

