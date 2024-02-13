function timespec(filename, maxiter, type, param)
% function timespec(filename, maxiter, type, param)

%
% $Id: timespec.m,v 1.7 2011/03/26 15:36:05 patrick Exp $
%
% Copyright (c) 2001-2011 Patrick Guio <patrick.guio@gmail.com>
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


t=[0:maxiter];

c=zeros(size(t)); 

switch type

	case 'sine'
		tau=param{1};
		c=50*(sin(2*pi*t/tau)+1);

	case 'cosine'
		tau=param{1};
		c=50*(cos(2*pi*t/tau)+1);

	case 'square'
		tau=param{1};
		ii=find(sin(2*pi*t/tau)>=0); 
		c(ii)=1000;
	
	case 'relaxsquare'
		tau=param{1};
		relax=param{2};
		ii=find(sin(2*pi*t/tau)>=0);
		c(ii)=1000;
		t0=t(1);
		for i=1:length(t)-1,
			if c(i)==1000 & c(i)==c(i+1),
				c(i)=c(i)*(1-exp(-(t(i)-t0)/relax));
			else
				t0=t(i);
			end
		end
		ii=find(sin(2*pi*t/tau)<0);
		for i=1:length(t)-1,
			if c(i)==0 & c(i)==c(i+1),
				c(i)=1000*(exp(-(t(i)-t0)/relax));
			else
				t0=t(i);
			end
		end

	case 'pulse'
		tmin=param{1};
		tmax=param{2};
		if length(param)==2,
			relax=0;
		else
			relax=param{3};
		end
		ii=find(t>=tmin & t<=tmax);
		c(ii)=1000;
		if (relax~=0)
			c(ii)=c(ii).*(1-exp(-(t(ii)-t(ii(1)))/relax));
			ii=find(t>tmax);
			c(ii)=1000.*(exp(-(t(ii)-t(ii(1)))/relax));
		end

end

c=round(c);

plot(t,c,'-o')
xlabel('iterations')
ylabel('coefficient')

fid=fopen(filename,'w');

fprintf(fid,'# maxiter=%d, type=''%s'', param={', maxiter, type);
for i=1:length(param)-1,
	fprintf(fid,'%d,', param{i});
end
fprintf(fid,'%d}\n',param{end});

fprintf(fid,'timespec=');
fprintf(fid,'%d,', c(1:end-1));
fprintf(fid,'%d', c(end));

fclose(fid);
