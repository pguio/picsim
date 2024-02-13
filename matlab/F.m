function F(m,s)
% function F(m,s)

%
% $Id: F.m,v 1.4 2011/03/26 15:36:04 patrick Exp $
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

x=linspace(0,10,100);
beta=1/sqrt(2*pi)*s*exp(-m^2/2/s^2)+m/2*(1+erf(m/sqrt(2)/s));

y=1/beta*x/sqrt(2*pi)/s.*exp(-(x-m).^2/2/s^2);

m0=d01gaf(x(:),y(:));
fprintf(1,'m0=%.4f\n', m0);

mu=1/beta*((m^2+s^2)/2*(1+erf(m/sqrt(2)/s))+m*s/sqrt(2*pi)*exp(-m.^2/2/s^2));
sigma=sqrt(1/beta*(m*(m^2+3*s^2)/2*(1+erf(m/sqrt(2)/s))+ ...
	s*(m^2+2*s^2)/sqrt(2*pi)*exp(-m.^2/2/s^2))-mu^2);

m1=d01gaf(x(:),x(:).*y(:));
m2=d01gaf(x(:),(x(:)-m1).^2.*y(:));

fprintf(1,'mu=%.4f (%.4f) sigma=%.4f (%.4f)\n', m1, mu,  sqrt(m2), sigma);

nb=500000;
r=randflux(m,s,1,nb);
hist(r,20);
fprintf(1,'mu=%.4f sigma=%.4f\n', mean(r), std(r));

