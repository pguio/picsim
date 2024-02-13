function field_plot(field,fieldname,scheduler,pos)
% function field_plot(field,fieldname,scheduler,pos)

%
% $Id: field_plot.m,v 1.5 2011/03/26 15:36:05 patrick Exp $
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

subplot(pos), 
%cla;
%set(gca,'box','on');
switch scheduler.dimr,
	case 2,
		imagesc(field'); 
		axis xy;
		colorbar('v');
		title(fieldname);
	case 3,
		imagesc(field(:,:,fix(end/2)+1)'); 
		axis xy;
		colorbar('v');
		title(fieldname);
end
drawnow
