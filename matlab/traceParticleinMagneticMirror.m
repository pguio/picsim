function traceParticleTrajectoryinMagneticDipole
% function traceParticleTrajectoryinMagneticDipole

%
% $Id: traceParticleinMagneticMirror.m,v 1.3 2015/11/26 17:00:16 patrick Exp $
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

tmax = 200;

options = odeset('RelTol',1e-10,'Stats','on');

%timespan
tspan = [0,tmax];
tspan = linspace(0,tmax,200);
Xo = [0; 0.5; 0; 0.1; 0.01; 0];
lbl1 = num2str(Xo');
lbl1 = sprintf('X=(%.2f,%.2f,%.2f) V=(%.2f,%.2f,%.2f)', Xo);
[t1,X1] = ode45(@dynamic,tspan,Xo,options);
Xo = [0; -0.3; 0; -.1; 0.02; 0];
lbl2 = num2str(Xo');
lbl2 = sprintf('X=(%.2f,%.2f,%.2f) V=(%.2f,%.2f,%.2f)', Xo);
[t2,X2] = ode15s(@dynamic,tspan,Xo,options);

clf
subplot(211)
plot(X1(:,1),sqrt(X1(:,2).^2+X1(:,3).^2), '-x', ...
     X2(:,1),sqrt(X2(:,2).^2+X2(:,3).^2), '-x')
xlabel('$x_\|$','interpreter','latex'),
ylabel('$x_\perp$','interpreter','latex'),
subplot(212)
plot(t1,sqrt(sum(X1(:,[4:6])'.^2)), t2, sqrt(sum(X2(:,[4:6])'.^2)))
xlabel('$t$','interpreter','latex'),
ylabel('$v$','interpreter','latex'),
set(gca,'xlim',[tspan(1), tspan(end)]);
set(gca,'ylim',[0.09 0.15]);
legend(lbl1,lbl2)


function [dr_dt]= dynamic(t,r)

[x,y,z,vx,vy,vz] = deal(r(1),r(2),r(3),r(4),r(5),r(6));

m = 1;
t1=pi/2;
RM1 = [3 0 0];
t2=pi/2;
RM2 = [-3 0 0];

M1 = [m*sin(t1), m*cos(t1), 0];
M2 = [m*sin(t2), m*cos(t2), 0];

B1 = dipoleMagneticField(M1,RM1,{x,y,z});
B2 = dipoleMagneticField(M2,RM2,{x,y,z});

[Bx,By,Bz] = deal(B1{1}+B2{1}, B1{2}+B2{2}, B1{3}+B2{3});

%[Bx, By, Bz]

dx = vx;
dy = vy;
dz = vz;

qOverm = 2; 

dvx = qOverm*(vy*Bz - vz*By);
dvy = qOverm*(vz*Bx - vx*Bz);
dvz = qOverm*(vx*By - vy*Bx);

dr_dt = [dx; dy; dz; dvx; dvy; dvz];



