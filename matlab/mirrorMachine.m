function B = mirrorMachine(R)
% function B = mirrorMachine(R)

%
% $Id: mirrorMachine.m,v 1.3 2016/06/06 14:33:59 patrick Exp $
%
% Copyright (c) 2010-2016 Patrick Guio <patrick.guio@gmail.com>
% All Rights Reserved.
%
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 3 of the License, or (at your
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.

m = 1;

t1=pi/2;
RM1 = [3 0];
t2=pi/2;
RM2 = [-3 0];
 
M1 = [m*sin(t1), m*cos(t1)];
M2 = [m*sin(t2), m*cos(t2)];

B1 = dipoleMagneticField(M1,RM1,R);
B2 = dipoleMagneticField(M2,RM2,R);

% B = { B1{1}+B2{1}, B1{2}+B2{2} };

Bx = B1{1}+B2{1};
By = B1{2}+B2{2};
Bz = zeros(size(B1{1}));

B = [Bx; By; Bz];


if 1,
Bn = sqrt(sum(B.^2));
if 1,
Bnmax  = 10;
ii = find(Bn > Bnmax);
B(:,ii) = B(:,ii)./(ones(3,1)*Bn(ii))*Bnmax;
Bn = sqrt(sum(B.^2));
end
else
Bn = sqrt(sum(B.^2));
B = B./(ones(3,1)*Bn);
end

X = R{1};
Y = R{2};

Bn = sqrt(sum(B.^2));
Bxu = B(1,:)./Bn;
Byu = B(2,:)./Bn;
clf
quiver(X,Y,Bxu,Byu,0.5)

fprintf(1,'min(B) %g max(B) %g\n', min(Bn), max(Bn));
pause(1)

