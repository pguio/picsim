function mirror
% function mirror

%
% $Id: mirror.m,v 1.1 2016/06/06 14:33:40 patrick Exp $
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

xlim = [-10 10];
ylim = [-3 3];
clim = [0 0.5];

m = 1;

x = linspace(-10,10,401);
y = linspace(-3,3,201);
[X,Y] = meshgrid(x,y);

% \vec{B}(\vec{m},\vec{r} ) = \frac{\mu_0}{4\pi r^3}
%    \left(3(\vec{m}\cdot\vec{\hat{r}}) \vec{\hat{r}} -\vec{m} \right) +
%    \frac{2\mu_0}{3} \vec{m} \delta^3(\vec{r})
%     


if 1

t1=pi/2;
RM1 = [5 0];
t2=pi/2;
RM2 = [-5 0];

else

t1=0;
RM1 = [3 0];
t2=pi;
RM2 = [-3 0];

end

M1 = [m*sin(t1), m*cos(t1)];
M2 = [m*sin(t2), m*cos(t2)];

B1 = dipoleMagneticField(M1,RM1,{X,Y});
B2 = dipoleMagneticField(M2,RM2,{X,Y});
B = { B1{1}+B2{1}, B1{2}+B2{2} };

alpha1 = dipoleMagneticPotential(M1,RM1,{X,Y});
alpha2 = dipoleMagneticPotential(M2,RM2,{X,Y});
alpha = alpha1+alpha2;


subplot(211), 
[c,h] = contour(X,Y,log10(sqrt(B{1}.^2+B{2}.^2)),40); 
%clabel(c,h)
axis equal
title('B isocontours')


subplot(212), 
[c,h] = contour(X,Y,log10(alpha),40);
%clabel(c,h)
axis equal
title('\alpha isocontours')
return

pause

subplot(211)

B = sqrt(Bx.^2+By.^2);
Bx = Bx./B;
By = By./B;

if 0
pcolor(X,Y,B)
shading flat
%axis equal
colorbar
set(gca,'xlim',xlim);
set(gca,'ylim',ylim);
set(gca,'clim',clim);
else
sy=linspace(-2,2,20);
sx=2*ones(size(sy));
streamline(stream2(X,Y,Bx,By,sx,sy))
sx=linspace(-10,10,20);
sy=zeros(size(sx));
streamline(stream2(X,Y,Bx,By,sx,sy))
%quiver(X,Y,Bx,By)
end

subplot(212) 

if 0
B = m./R.^3.*sqrt(1+3*X.^2./R.^2);
else
B = m./R.^3.*sqrt(1+3*Y.^2./R.^2);
end

pcolor(X,Y,B)
shading flat
%axis equal
colorbar
set(gca,'xlim',xlim);
set(gca,'ylim',ylim);
set(gca,'clim',clim);


function B = magdip(m,rm, r)

if length(m) ~= length(r)
  error('length(m) should be identical to length(r)');
end

[X, Y] = deal(r{:});
X = X-rm(1);
Y = Y-rm(2);
R = sqrt(X.^2 + Y.^2);

mDotRu = (m(1)*X + m(2)*Y)./R;

B{1} = 1./R.^3.*(3*mDotRu.*X./R - m(1));
B{2} = 1./R.^3.*(3*mDotRu.*Y./R - m(2));


function a = magpot(m,rm,r)

if length(m) ~= length(r)
  error('length(m) should be identical to length(r)');
end

[X, Y] = deal(r{:});

X = X-rm(1);
Y = Y-rm(2);

R2 = X.^2 + Y.^2;
R = sqrt(X.^2 + Y.^2);

mCrossRu2 = (m(1)*Y - m(2)*X).^2./R2;

a = mCrossRu2./R;

size(a)
size(a>0)
size(isfinite(a))

function a = equipot(m, rm, r, a0)

if length(m) ~= length(r)
  error('length(m) should be identical to length(r)');
end

[X, Y] = deal(r{:});
%X = X-rm(1);
%Y = Y-rm(2);
R2 = X.^2 + Y.^2;
R = sqrt(X.^2 + Y.^2);

mCrossRu2 = (m(1)*Y - m(2)*X).^2./R2;

a{1} = a0*mCrossRu2.*X./R+rm(1);
a{2} = a0*mCrossRu2.*Y./R+rm(2);

