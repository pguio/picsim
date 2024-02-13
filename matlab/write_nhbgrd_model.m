function write_nhbgrd_model(nhbgrd,filename)
% function write_nhbgrd_model(nhbgrd,filename)

%
% $Id: write_nhbgrd_model.m,v 1.3 2017/12/23 18:46:33 patrick Exp $
%
% Copyright (c) 2009-2011 Patrick Guio <patrick.guio@gmail.com>
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


switch ndims(nhbgrd),

  case 2,
	  [nx,ny] = size(nhbgrd);
	  filename = [filename '.nh2d'];
		fid = fopen(filename, 'w');
		fprintf(fid,'%d %d\n', nx, ny);
		% save colum by column, i.e. first coordinate change fastest
		for j=1:ny,
		  fprintf(fid,'%.10g ', nhbgrd(:,j));
			fprintf(fid,'\n');
		end
		fclose(fid);

	case 3,
	  [nx,ny,nz] = size(nhbgrd);
	  filename = [filename '.nh3d'];
		fid = fopen(filename, 'w');
		fprintf(fid,'%d %d %d\n', nx, ny, nz);
		for k=1:nz,
		  for j=1:ny,
		    fprintf(fid,'%10.8f ', nhbgrd(:,j,k));
			  fprintf(fid,'\n');
			end
		end
		fclose(fid);

end

fprintf(1,'Testing reading... ');
nhbgrd1 = read_nhbgrd_model(filename);
fprintf(1,'ok\n');

diff = abs(nhbgrd1-nhbgrd)./nhbgrd;
fprintf(1,'Relative error min = %.2g %% / max = %.2g %% \n', ...
        100*[min(diff(:)), max(diff(:))]);


