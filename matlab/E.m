function E(lambda)
% function E(lambda)

%
% $Id: E.m,v 1.4 2011/03/26 15:36:04 patrick Exp $
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

x=linspace(0,30,100);
y=1/lambda*exp(-x/lambda);

m0=d01gaf(x(:),y(:));
fprintf(1,'m0=%.4f\n', m0);

mu=lambda;
sigma=lambda;

m1=d01gaf(x(:),x(:).*y(:));
m2=d01gaf(x(:),(x(:)-m1).^2.*y(:));

fprintf(1,'mu=%.4f (%.4f) sigma=%.4f (%.4f)\n', m1, mu,  sqrt(m2), sigma);

nb=500000;
r=lambda*randexp(1,nb);
hist(r,20);
fprintf(1,'mu=%.4f sigma=%.4f\n', mean(r), std(r));

