function [scheduler,species]=picsim_init2_2d(varargin)
% function [scheduler,species]=picsim_init2_2d(varargin)

%
% $Id: picsim_init2_2d.m,v 1.5 2011/03/26 15:36:05 patrick Exp $
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
scheduler = scheduler_par_parsing(scheduler,varargin{:});
% scheduler initialisation
scheduler = scheduler_init(scheduler);

% beam definition file
beam = beam_def(scheduler);
% beam parameter parsing
beam = beam_par_parsing(beam,varargin{:});
% beam initialisation
beam = beam_init(beam,scheduler);

species = { beam };
