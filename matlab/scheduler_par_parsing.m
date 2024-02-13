function scheduler=scheduler_par_parsing(scheduler,varargin)
% function scheduler=scheduler_par_parsing(scheduler,varargin)

%
% $Id: scheduler_par_parsing.m,v 1.5 2011/03/26 15:36:05 patrick Exp $
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

if isempty(varargin),
	return
end

argc=length(varargin);

if rem(argc,2)~=0,
  error('Extra arguments should be in pair parameter/value')
end

for i=1:2:argc,

	param=varargin{i};
	val=varargin{i+1};
	switch lower(param)
		case 'dimr', scheduler.dimr=val;
		case 'dimv', scheduler.dimv=val;

		case 'ustate', scheduler.ustate=val;
		case 'nstate', scheduler.nstate=val;

		case 'dt', scheduler.dt=val;
		case 'maxiter', scheduler.maxiter=val;

		case 'mass', scheduler.mass=val;
		case 'lambda', scheduler.Lambda=val;
		case 'temp', scheduler.temp=val;

		case 'nx', scheduler.nx=val;
		case 'ny', scheduler.ny=val;
		case 'nz', scheduler.nz=val;
		case 'xa', scheduler.rmin(1)=val;
		case 'xb', scheduler.rmax(1)=val;
		case 'yc', scheduler.rmin(2)=val;
		case 'yd', scheduler.rmax(2)=val;
		case 'ze', scheduler.rmin(3)=val;
		case 'zf', scheduler.rmax(3)=val;

		case 'electrostatic', scheduler.electrostatic=val;
		case 'electromagnetic', scheduler.electromagnetic=val;
		case 'gravitation', scheduler.gravitation=val;

		case 'extmagneticfield', scheduler.extMagneticField = val;

		case 'tracer', scheduler.tracer = val;
	end

end
