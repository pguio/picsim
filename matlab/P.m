function P(m,s,alpha)
% function P(m,s,alpha)

%
% $Id: P.m,v 1.7 2011/03/26 15:36:05 patrick Exp $
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

x=linspace(0,1,100);

%if alpha>=0
	k=2*alpha/(erf((1-m)/sqrt(2)/s)+erf(m/sqrt(2)/s));
	beta=1+k/2*(erf((1-m)/sqrt(2)/s)+erf(m/sqrt(2)/s));
	y=1/beta*(1+k/sqrt(2*pi)/s*exp(-(x-m).^2/2/s^2));
%else
%	k=2*alpha/(2*alpha/sqrt(2*pi)/s+erf((1-m)/sqrt(2)/s)-erf(-m/sqrt(2)/s));
%	beta=1-k/sqrt(2*pi)/s+k/2*(erf((1-m)/sqrt(2)/s)+erf(m/sqrt(2)/s));
%	y=1/beta*(1+k/sqrt(2*pi)/s*(exp(-(x-m).^2/2/s^2)-1));
%end

m0=d01gaf(x(:),y(:));
fprintf(1,'m0=%.4f\n', m0);

%if alpha>=0
	mu = 1/beta*(1/2-k*s/sqrt(2*pi)*(exp(-(1-m)^2/2/s^2)-exp(-m^2/2/s^2))+ ...
	    k*m/2*(erf((1-m)/sqrt(2)/s)+erf(m/sqrt(2)/s)));
	sigma = sqrt(1/beta*(1/3- ...
	    k/sqrt(2*pi)*s*((1+m)*exp(-(1-m)^2/2/s^2)-m*exp(-m^2/2/s^2))+ ...
			k*(m^2+s^2)/2*(erf((1-m)/sqrt(2)/s)+erf(m/sqrt(2)/s)))-mu^2);
%else
%	mu = 1/beta*((1-k/sqrt(2*pi)/s)/2- ...
%			k*s/sqrt(2*pi)*(exp(-(1-m)^2/2/s^2)-exp(-m^2/2/s^2))+ ...
%			k*m/2*(erf((1-m)/sqrt(2)/s)+erf(m/sqrt(2)/s)));
%	sigma = sqrt(1/beta*((1-k/sqrt(2*pi)/s)/3- ...
%	    k/sqrt(2*pi)*s*((1+m)*exp(-(1-m)^2/2/s^2)-m*exp(-m^2/2/s^2))+ ...
%			k*(m^2+s^2)/2*(erf((1-m)/sqrt(2)/s)+erf(m/sqrt(2)/s)))-mu^2);
%end

m1=d01gaf(x(:),x(:).*y(:));
m2=d01gaf(x(:),(x(:)-m1).^2.*y(:));

fprintf(1,'mu=%.4f (%.4f) sigma=%.4f (%.4f)\n', m1, mu,  sqrt(m2), sigma);

nb=500000;
r=randperturb(m,s,alpha,1,nb);
hist(r,20);
fprintf(1,'mu=%.4f sigma=%.4f\n', mean(r), std(r));

