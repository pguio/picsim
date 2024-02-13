function [scheduler,species]=picsim_init4(varargin)
% function [scheduler,species]=picsim_init4(varargin)

%
% $Id: picsim_init4_2d.m,v 1.4 2011/03/26 15:36:05 patrick Exp $
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

% default scheduler definition setup
scheduler = scheduler2d_def;
% scheduler parameter parsing
if ~isempty(varargin),
  scheduler = scheduler_par_parsing(scheduler,varargin{:});
else
  opts = {'maxiter',500, ...
	        'lambda',20, ...
          'nx', 113, 'ny', 81, ...
	        'xa',-5, 'xb', 5, ...
          'yc', -2, 'yd', 2, ...
					'electrostatic',0, ...
					'electromagnetic',2, ...
					'extMagneticField',@mirrorMachine, ...
					'tracer', true};
  scheduler = scheduler_par_parsing(scheduler, opts{:});
end
% scheduler initialisation
scheduler = scheduler_init(scheduler);

% background definition file
background = backgrd_def(scheduler);
% background parameter parsing
if ~isempty(varargin),
  background = backgrd_par_parsing(background,varargin{:});
else
  opts = {'lambda',20, ....
	        'tempx',.01,'tempy',.01,'tempz',.01, ...
	        'bc', 0 };
	background = backgrd_par_parsing(background, opts{:});
end
% background initialisation
background = backgrd_init(background,scheduler);

species = { background };

