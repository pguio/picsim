function rotation
% function rotation

%
% $Id: rotation.m,v 1.3 2011/03/26 15:36:05 patrick Exp $
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

B=[0 1 0];
dt=0.2;
Z=1;
m=1;

T=Z*B/m*dt/2;
S=2*T/(1+sum(T.*T));

v1=[1 0 0];
fprintf(1,'|v1|=%f\n', norm(v1));

v=v1+cross(v1,T);
fprintf(1,'|v| =%f\n', norm(v));

v2=v1+cross(v,S);
fprintf(1,'|v2|=%f\n', norm(v2));


fprintf(1,'acos(v1,v2)=%f\n', acos(sum(v1.*v2)/(norm(v1)*norm(v2))));
fprintf(1,'asin(v1,v2)=%f\n', asin(norm(cross(v1,v2))/(norm(v1)*norm(v2))));
fprintf(1,'atan(v1,v2)=%f\n', 2*atan(norm(T)));

clf
h=line([0 v1(1)], [0 v1(3)]); 
set(h,'linestyle','-');
h=line([0 v(1)], [0 v(3)]); 
set(h,'linestyle','--');
h=line([0 v2(1)], [0 v2(3)]); 
set(h,'linestyle','-.');
axis equal


v1=[1 0 0];
fprintf(1,'|v1|=%f\n', norm(v1));

v2=v1+cross(v1,T);
v2=v2*norm(v1)/norm(v2);
fprintf(1,'|v| =%f\n', norm(v2));

fprintf(1,'acos(v1,v2)=%f\n', acos(sum(v1.*v2)/(norm(v1)*norm(v2))));
fprintf(1,'asin(v1,v2)=%f\n', asin(norm(cross(v1,v2))/(norm(v1)*norm(v2))));

fprintf(1,'2*acos(v1,v2)=%f\n', 2*acos(sum(v1.*v2)/(norm(v1)*norm(v2))));
fprintf(1,'2*asin(v1,v2)=%f\n', 2*asin(norm(cross(v1,v2))/(norm(v1)*norm(v2))));

h=line([0 v(1)], [0 v(3)]);
set(h,'linestyle','--');
h=line([0 v2(1)], [0 v2(3)]);
set(h,'linestyle','-.');
axis equal

