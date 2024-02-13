function field=readhdf(filename,fieldname,subdim,seeds)
% function field=readhdf(filename,fieldname,subdim,seeds)

%
% $Id: readhdf.m,v 1.14 2011/03/26 15:36:05 patrick Exp $
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

if ~exist('subdim','var')
	subdim=[];
end

if ~exist('seeds','var') | isempty(seeds),

	field=readone(filename,fieldname,subdim);

else

	file = [filename num2str(seeds(1))];
	field=readone(file,fieldname,subdim);

	var=zeros([length(seeds) size(field.var)]);

	for i=1:length(seeds),
		file = [filename num2str(seeds(i))];
		fprintf(1,'Loading %s\n',file);
		f=readone(file,fieldname,subdim);
		var(i,:)=f.var(:);
	end

	field.var=reshape(mean(var,1),size(field.var));
	field.var_std=reshape(std(var,1),size(field.var));

end


function field=readone(filename,fieldname,subdim)

if isempty(findstr(filename,'.hdf')),
	filename=[filename '.hdf'];
end

sd_id=hdfsd('start',filename,'rdonly');
idx = hdfsd('nametoindex',sd_id',fieldname);
sds_id = hdfsd('select',sd_id,idx);
[name,rank,dimsizes,data_type,nattrs,status] = hdfsd('getinfo',sds_id);
fprintf(1,'Variable `%s'' has %d dimensions\n', name, rank);

if ~exist('subdim') | isempty(subdim),
	start=zeros(size(dimsizes));
	stride=ones(size(dimsizes));
	edge=dimsizes;
else
	for i=1:length(subdim),
		if isempty(subdim{i})
			start(i)=0;
			stride(i)=1;
			edge(i)=dimsizes(i);
		else
			start(i)=subdim{i}(1)-1;
			if length(subdim{i})>1,
				stride(i)=subdim{i}(2)-subdim{i}(1);
			else
				stride(i)=1;
			end
			edge(i)=length(subdim{i});
		end
	end
end

%fprintf(1,'start=%d stride=%d edge=%d\n',[start;stride;edge])

[var,status]=hdfsd('readdata',sds_id,start,stride,edge);
if rank>1,
	field.var=double(permute(var,rank:-1:1));
else
	field.var=double(var);
end
field.varname=convert(name);

for dim_number=0:rank-1,
	dim_id = hdfsd('getdimid',sds_id,dim_number);
end

for dim_number=0:rank-1,
	dim_id = hdfsd('getdimid',sds_id,dim_number);
	[vardim,status]=hdfsd('getdimscale',dim_id);
	field.dims{dim_number+1}=double(vardim);
	if exist('subdim') & ~isempty(subdim) & ~isempty(subdim{dim_number+1}),
		field.dims{dim_number+1}=field.dims{dim_number+1}(subdim{dim_number+1});
	end

	[name,count,data_type,nattrs,status] = hdfsd('diminfo',dim_id);
	fprintf(1,'--> %4.d elements in dimension `%s''\n', count,name);
	field.dimsname{dim_number+1}=convert(name);
end

hdfsd('endaccess',sds_id);
hdfsd('end',sd_id);

tt=convert(fieldname);

function title=convert(tt)

i=findstr(tt,'-');
if ~isempty(i)
	tt=tt(1:i-1);
end

switch(tt)
case 'potential', title='\phi';
case 'density', title='\rho';
case 'Ex-field', title='E_x';
case 'Ey-field', title='E_y';
case 'x-velocity', title='u_x'; 
case 'y-velocity', title='u_y'; 
case 'z-velocity', title='u_z'; 
case 'x-temperature', title='T_x'; 
case 'y-temperature', title='T_y'; 
case 'z-temperature', title='T_z'; 
case 'x-vx', title='Phase space x-v_x'; 
case 'x-vy', title='Phase space x-v_y'; 
case 'x-vz', title='Phase space x-v_z'; 
case 'y-vx', title='Phase space y-v_x'; 
case 'y-vy', title='Phase space y-v_y'; 
case 'y-vz', title='Phase space y-v_z'; 
case 'ese', title='ESE';
case 'time', title='\tau';
case 'Time-spectra', title='\tau';
case 'frequency', title='frq';
case 'Time', title='\tau';
case 'avSkw2', title='<|N_e|^2>';
case 'K', title='k';
otherwise, title=tt; %error('unknown title');
end
