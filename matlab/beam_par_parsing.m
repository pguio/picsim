function particles=beam_par_parsing(particles,varargin)
% function particles=beam_par_parsing(particles,varargin)

%
% $Id: beam_par_parsing.m,v 1.4 2011/03/26 15:36:05 patrick Exp $
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
		case 'charge', particles.charge=val;
		case 'mass', particles.mass=val;
		case 'lambda', particles.Lambda=val;
		case 'lambda', particles.Lambda=val;

		case 'beamdir', particles.beamdir=val;
		case 'u', particles.u(:)=val;

		case 'tempx', particles.temp(1)=val;
		case 'tempy', particles.temp(2)=val;
		case 'tempz', particles.temp(3)=val;

	end

end
