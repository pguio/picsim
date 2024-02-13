function nhbgrd = read_nhbgrd_model(filename)
% function nhbgrd = read_nhbgrd_model(filename)

%
% $Id: read_nhbgrd_model.m,v 1.2 2011/03/26 15:36:05 patrick Exp $
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



if ~isempty(strfind(filename,'2d'))
  dim = 2;
elseif ~isempty(strfind(filename,'3d'))
  dim = 3;
end

switch dim,

  case 2,
		fid = fopen(filename, 'r');
		[a,count] = fscanf(fid,'%d', 2);
		nx = a(1);
		ny = a(2);
		nhbgrd = zeros(nx,ny);
		for j=1:ny,
		  [a,count] = fscanf(fid,'%f', nx);
			nhbgrd(:,j) = a;
		end
		%fclose(fid);

	case 3,
		fid = fopen(filename, 'r');
		[a,count] = fscanf(fid,'%d %d\n', 1);
		nx = a(1);
		ny = a(2);
		nz = a(3);
		nhbgrd = zeros(nx,ny,nz);
		for k=1:nz,
		  for j=1:ny,
			  [a,count] = fscanf(fid,'%10.8f ', nx);
				nhbgrd(:,j,k) = a;
			  [a,count] = fscanf(fid,'\n', 1);
		end
			end
		end
		%fclose(fid);

end
