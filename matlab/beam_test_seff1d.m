function beam_test_seff(xmin,xmax,n,m,s,dx,seed)
% function beam_test_seff(xmin,xmax,n,m,s,dx,seed)

%
% $Id: beam_test_seff1d.m,v 1.11 2011/03/26 15:36:05 patrick Exp $
%
% Copyright (c) 2004-2011 Patrick Guio <patrick.guio@gmail.com>
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

x=xmin:dx:xmax;
xp=x+0.5*dx;
xm=x-0.5*dx;

% 
ipmx=max(find(m>=xp));
imx=max(find(m>=x));
fprintf(1, 'xp(%d)=%+g xp(%d)=%+g\n', [ipmx xp(ipmx) ipmx+1 xp(ipmx+1)]);
fprintf(1, 'x(%d) =%+g x(%d) =%+g\n', [imx x(imx) imx+1 x(imx+1)]);

% cumulated prob [xmin xmax]
beta0=igauss(x(end),m,s)-igauss(x(1),m,s);
fprintf(1, 'beta0=%g\n', beta0);

% n is number in the interval around m
%N=fix(n/(igauss(xp(ipmx+1),m,s)-igauss(xp(ipmx),m,s))*beta0);
N=fix(n/max(igauss(xp,m,s)-igauss(xm,m,s))*beta0);
fprintf('n=%d N=%d N/n=%g\n', n, N, N/n);

% Generate truncated random gaussian(m,s) between [xmin xmax[
randn('seed',seed);
Y=s*randn(1,N)+m;
while ~isempty(find(Y<x(1) | Y>=x(end))),
	ii=find(Y<x(1) | Y>=x(end));
	Y(ii) = s*randn(1,length(ii))+m;
end

% histogram with bins centered on x
% zeroth-order field weighting
Nx1=zeros(size(x));
for i=1:length(x),
  ii=find(Y>=xm(i) & Y<xp(i));
  Nx1(i) = length(ii);
end

% constant approximation with bins centered on x
Nx2=gauss(x,m,s)*dx;
beta1=sum(Nx2);
fprintf(1, 'beta1=%g\n', beta1);
fprintf('n=%d N=%d N/n=%g\n', n, fix(n/max(Nx2)*beta1), beta1/max(Nx2));
Nx2=Nx2*n/max(Nx2);

% integrated probability with bins centered on x
Nx3=N*(igauss(xp,m,s)-igauss(xm,m,s))/beta0;

%plot(x, Nx1, '-o', x, Nx2, x, Nx3);
%return

% linear interpolation
% first-order field weighting
Nx4=zeros(size(x));
for i=1:length(x)-1,
	ii=find(Y>=x(i) & Y<x(i+1));
	Nx4(i) = Nx4(i)+sum(1-(Y(ii)-x(i))/dx);
	Nx4(i+1)= Nx4(i+1)+sum(1-(x(i+1)-Y(ii))/dx);
end

%plot(x, Nx1, '-o', x, Nx4, '-x');
%return

Nx5=zeros(size(x));

Nx5(1)=fp(x(2),x(1),m,s,dx)-fp(x(1),x(1),m,s,dx);
Nx5(2:end-1)=fp(x(3:end),x(2:end-1),m,s,dx)-fp(x(2:end-1),x(2:end-1),m,s,dx)+...
fm(x(2:end-1),x(2:end-1),m,s,dx)-fm(x(1:end-2),x(2:end-1),m,s,dx);
Nx5(end)=fm(x(end),x(end),m,s,dx)-fm(x(end-1),x(end),m,s,dx);

beta2=sum(Nx5);
fprintf(1, 'beta2=%g\n', beta2);
fprintf('n=%d N=%d N/n=%g\n', n, fix(n/max(Nx5)*beta2), beta2/max(Nx5));
Nx5=Nx5*n/max(Nx5);

plot(x, Nx1, 'o', x, Nx2, x, Nx3, x, Nx4, 'x', x, Nx5,'markersize',5);

fprintf(1,'%d %.0f %.0f %.0f %.0f\n', ...
	max(Nx1), max(Nx2), max(Nx3), max(Nx4), max(Nx5));


function y=gauss(x,m,s)
y=1/sqrt(2*pi)/s*exp(-(x-m).^2/2/s^2);


function y=igauss(x,m,s)
% y=\int gauss(x,m,s)
y=1/2*erf((x-m)/sqrt(2)/s);


function y=fm(x,x0,m,s,dx)
% y=\int 1/dx*( (x-x0)/dx+1)/sqrt(2*pi)/s.*exp(-(x-m).^2/2/s^2);
y=1/dx^2*(-s/sqrt(2*pi)*exp(-(x-m).^2/2/s^2)+1/2*(m-x0+dx).*erf((x-m)/sqrt(2)/s));
%y=1/2*(-sqrt(2)*s*exp(-(x-m).^2/2/s^2)+m*sqrt(pi)*erf((x-m)/sqrt(2)/s)-x0*sqrt(pi).*erf((x-m)/sqrt(2)/s)+sqrt(pi)*erf((x-m)/sqrt(2)/s)*dx)/dx^2/sqrt(pi);

function y=fp(x,x0,m,s,dx)
% y=\int 1/dx*(-(x-x0)/dx+1)/sqrt(2*pi)/s.*exp(-(x-m).^2/2/s^2);
y=1/dx^2*(s/sqrt(2*pi)*exp(-(x-m).^2/2/s^2)-1/2*(m-x0-dx).*erf((x-m)/sqrt(2)/s));
%y=1/2*(sqrt(2)*s*exp(-(x-m).^2/2/s^2)-m*sqrt(pi)*erf((x-m)/sqrt(2)/s)+x0*sqrt(pi).*erf((x-m)/sqrt(2)/s)+sqrt(pi)*erf((x-m)/sqrt(2)/s)*dx)/dx^2/sqrt(pi);


