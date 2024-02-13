function traceParticleTrajectoryinMagneticDipole
% function traceParticleTrajectoryinMagneticDipole

%
% $Id: traceParticleinMagneticDipole.m,v 1.4 2015/11/26 17:00:16 patrick Exp $
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

tmax = 100;
options = odeset('AbsTol',1e-14,'RelTol',1e-12,'Stats','on','vectorized','on');
tspan = [0,tmax];


for i=1:4,
tspans = linspace(0,tmax,1000);
Xo = [1; 0; 0; .01; 0.2; 0];
X=[Xo,Xo,Xo];
%size(X)
%dX=dynamic3d([0,0,0],X)
%return
switch i,
case 1, [t1,X1] = ode45(@dynamic3d,tspans,Xo,options); solver='ode45';
case 2, [t1,X1] = ode15s(@dynamic3d,tspans,Xo,options); solver='ode15s';
case 3, [t1,X1] = ode87(@dynamic3d,tspan,Xo,options); solver='ode87';
case 4, [t1,X1] = ode78(@dynamic3d,tspan,Xo,0,options.RelTol,0); solver='ode78';
end
size(t1), size(X1)
%size(X1)
clf
E = X1(:,4).^2+X1(:,5).^2+X1(:,6).^2;
meanE = mean(E);
stdE = std(E);
fprintf(1,'mean E=%4g, sigma E=%4g\n', meanE,stdE);
subplot(211), plot(sqrt(X1(:,1).^2+X1(:,3).^2),X1(:,2)), 
xlabel('$\sqrt{x^2+z^2}$','interpreter','latex'),
ylabel('$y$','interpreter','latex'),
title(solver)
legend({sprintf('$X=(%.0f,%.0f,%.0f), V=(%.2f,%.2f,%.2f)$', Xo)},...
      'interpreter','latex');
set(gca,'xlim',[0 2])
set(gca,'ylim',[-0.5 0.5])
%axis equal,
subplot(212), plot(X1(:,1),X1(:,3)), 
xlabel('$x$','interpreter','latex'),
ylabel('$z$','interpreter','latex'),
title(sprintf('$m_E=$%4g, $\\sigma_E=$%4g', meanE, stdE),'interpreter','latex');
set(gca,'xlim',[-1.5 1.5])
set(gca,'ylim',[-1.5 1.5])
%axis equal
%subplot(313), plot(t1,E)
%xlabel('t')
%ylabel('|V|^2')
drawnow
orient tall
if i==1, print(gcf,'-dpsc2','cmp.ps');
else, print(gcf,'-dpsc2','-append','cmp.ps');
end
if 0
clf
plot3(X1(:,1),X1(:,2),X1(:,3))
set(gca,'xlim',[-1.5 1.5])
set(gca,'ylim',[-1.5 1.5])
set(gca,'zlim',[-1.5 1.5])
pause
end
end
%return

return
R1 = sqrt(X1(:,1).^2+X1(:,3).^2);
Z1 = X1(:,2);
E1 = 100*sum(X1(:,[4:6])'.^2)/sum(Xo(4:6).^2);

Xo = [1; 0; .01; 0.2; 0];
[t2,X2] = ode45(@dynamic2d,tspan,Xo,options);
%[t2,X2] = ode87(@dynamic2d,tspan,Xo,options);
x2 = X2(:,1);
y2 = X2(:,2);
E2 = 100*sum(X2(:,[3:5])'.^2)/sum(Xo(3:5).^2);

subplot(211)
plot(R1, Z1,'-x', x2, y2,'-x')
subplot(212)
plot(t1,E1,t2,E2)
set(gca,'xlim',tspan);
set(gca,'ylim',[80 120]);

return



Xo = [1; 0; 0; -.02; 0.2; 0];
[t2,X2] = ode45(@dynamic3d,tspan,Xo,options);
R2 = sqrt(X2(:,1).^2+X2(:,3).^2);
Z2 = X2(:,2);
E2 = 100*sum(X2(:,[4:6])'.^2)/sum(Xo(4:6).^2);

clf
subplot(211)
plot(R1, Z1,'-x', R2, Z2,'-x')
subplot(212)
plot(t1,E1,t2,E2)
set(gca,'xlim',tspan);
set(gca,'ylim',[80 120]);


function dr_dt = dynamic3d(t,r)

dr_dt = zeros(size(r));
[x,y,z,vx,vy,vz] = deal(r(1,:),r(2,:),r(3,:),r(4,:),r(5,:),r(6,:));

m = 1;
B = dipoleMagneticField([0, m, 0], [0, 0, 0], {x,y,z});
[Bx,By,Bz] = deal(B{:});

Z = 1; 

dvx = Z*(vy.*Bz - vz.*By);
dvy = Z*(vz.*Bx - vx.*Bz);
dvz = Z*(vx.*By - vy.*Bx);

dr_dt = [vx;vy;vz;dvx;dvy;dvz];

[r(1,:),r(2,:),r(3,:),r(4,:),r(5,:),r(6,:)] = deal(vx,vy,vz,dvx,dvy,dvz);
dr_dt = r;


function dr_dt = dynamic2d(t,r)

[x,y,vx,vy,vz] = deal(r(1,:),r(2,:),r(3,:),r(4,:),r(5,:));

m = 1;
B = dipoleMagneticField([0, m], [0, 0], {x,y});
[Bx,By] = deal(B{:});

Z = 1;

dvx = -Z*vz.*By;
dvy =  Z*vz.*Bx;
dvz =  Z*(vx.*By - vy.*Bx);

dr_dt = [vx;vy;dvx;dvy;dvz];

