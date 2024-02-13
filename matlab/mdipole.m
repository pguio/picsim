function mdipole
% function mdipole

%
% $Id: mdipole.m,v 1.1 2016/06/06 14:31:58 patrick Exp $
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

xlim = [-3 3];
ylim = [-3 3];
clim = [0 5];

close all

subplot(211)

% dipole moment [A m^2] or [J T^-1]
m = 1;

% distance from centre
r = linspace(0.2,3,100);
% magnetic colatitude \theta is 0 along the dipole's axis and
% 90 deg in the plane perpendicular to its axis
theta = linspace(0,pi,100);
[R,T] = meshgrid(r,theta);
% magnetic latitude \lambda = 90 deg − \theta  
lambda = pi/2-theta;
L = pi/2-T;
X = R.*cos(L);
Y = R.*sin(L);
B = m./R.^3.*sqrt(1+3*sin(L).^2);
pcolor(X,Y,B)
shading flat
%axis equal 
colorbar
set(gca,'xlim',xlim);
set(gca,'ylim',ylim);
set(gca,'clim',clim);
t=linspace(-pi/2,pi/2,100);
%re = linspace(1,5,5);
re = logspace(0,2,6);
xi = cell(size(re));
yi = cell(size(re));
hold on
for i=1:length(re),
  r = re(i)*sin(pi/2-t).^2;
  xi{i} = r.*cos(t);
  yi{i} = r.*sin(t);
	plot(xi{i},yi{i},'k')
end
hold off


subplot(212)

x = linspace(-3,3,300);
y = linspace(-3,3,300);
[X,Y] = meshgrid(x,y);
R = sqrt(X.^2+Y.^2);

% \vec{B}(\vec{m},\vec{r} ) = \frac{\mu_0}{4\pi r^3}
%    \left(3(\vec{m}\cdot\vec{\hat{r}}) \vec{\hat{r}} -\vec{m} \right) +
%    \frac{2\mu_0}{3} \vec{m} \delta^3(\vec{r})
%    
Bx = m./R.^3.*(3*Y./R.*X./R);
By = m./R.^3.*(3*Y./R.*Y./R-1);

tilt=pi/5;
M = [m*sin(tilt), m*cos(tilt)];
RM = [0 1];
Bv = dipoleMagneticField(M,RM,{X,Y});

[Bx,By] = deal(Bv{:});

alpha = dipoleMagneticPotential(M,RM,{X,Y});

if 0
subplot(221), contour(X,Y,log10(sqrt(Bx.^2+By.^2)),20); axis equal
subplot(222), contour(X,Y,log10(alpha)), axis equal
end


t=linspace(0,2*pi,200);
for i=1:length(xi),
  a{i} = dipoleMagneticEquiPotential(M, RM, {cos(t), sin(t) }, re(i));
end

B = sqrt(Bx.^2+By.^2);

%B = m./R.^3.*sqrt(1+3*Y.^2./R.^2);

pcolor(X,Y,B)
shading flat
axis equal
colorbar
set(gca,'xlim',xlim);
set(gca,'ylim',ylim);
set(gca,'clim',clim);

hold on
for i=1:length(xi),
  r = a{i};
	plot(r{1},r{2},'k')
end
hold off

figure

s = 6;

Xu = X(1:s:end,1:s:end);
Yu = Y(1:s:end,1:s:end);
Bxu = Bx(1:s:end,1:s:end)./B(1:s:end,1:s:end);
Byu = By(1:s:end,1:s:end)./B(1:s:end,1:s:end);


Ru = sqrt((Xu-RM(1)).^2+(Yu-RM(2)).^2);
ij = find(Ru>1);
[i,j] = ind2sub(size(Ru),ij);

Bx = Bx./B;
By = By./B;

if 1
quiver(Xu(ij),Yu(ij),Bxu(ij),Byu(ij),0.25)
set(gca,'xlim',xlim);
set(gca,'ylim',ylim);
end

x = X(1:s:end,1:s:end);
y = Y(1:s:end,1:s:end);
px = Bx(1:s:end,1:s:end);
py = By(1:s:end,1:s:end);

if 0
quiver(x,y,px,py,0.25)
set(gca,'xlim',xlim);
set(gca,'ylim',ylim);
end

hold on
for i=1:length(xi),
  r = a{i};
	plot(r{1},r{2},'k')
end
hold off

