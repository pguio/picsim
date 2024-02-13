function [scheduler,species]=picsim_init1_3d(varargin)
% function [scheduler,species]=picsim_init1_3d(varargin)

%
% $Id: picsim_init1_3d.m,v 1.5 2011/03/26 15:36:05 patrick Exp $
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
scheduler = scheduler3d_def;
% scheduler parameter parsing
scheduler = scheduler_par_parsing(scheduler,varargin{:});
% scheduler initialisation
scheduler = scheduler_init(scheduler);

% background definition file
background = backgrd_def(scheduler);
% background parameter parsing
background = backgrd_par_parsing(background,varargin{:});
% background initialisation
background = backgrd_init(background,scheduler);

species = { background };
