function U(a,b)
% function U(a,b)

%
% $Id: U.m,v 1.5 2011/03/26 15:36:05 patrick Exp $
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

x=linspace(a,b,100);
y=repmat(1/(b-a),size(x));

m0=d01gaf(x(:),y(:));
fprintf(1,'m0=%.4f\n', m0);

mu=(b+a)/2;
sigma=(b-a)/sqrt(12);

m1=d01gaf(x(:),x(:).*y(:));
m2=d01gaf(x(:),(x(:)-m1).^2.*y(:));

fprintf(1,'mu=%.4f (%.4f) sigma=%.4f (%.4f)\n', m1, mu,  sqrt(m2), sigma);

nb=500000;
r=(b-a)*rand(1,nb)+a;
hist(r,20);
fprintf(1,'mu=%.4f sigma=%.4f\n', mean(r), std(r));
