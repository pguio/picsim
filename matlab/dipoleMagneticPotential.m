function A = dipoleMagneticPotential(m, rm, r)
% function A = dipoleMagneticPotential(m, rm, r)
%
%    m  : magnetic moment vector
%    rm : position of the magnetic moment
%    r  : cartesion coordinates where to calculate magnetic field
%    A  : dipole magnetic potential (same size as r)

%
% $Id: dipoleMagneticPotential.m,v 1.2 2011/03/26 15:36:05 patrick Exp $
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

if length(m) ~= length(r)
  error('length(m) should be identical to length(r)');
end

% \vec{B}(\vec{m},\vec{r} ) = \frac{\mu_0}{4\pi r^3}
%    \left(3(\vec{m}\cdot\vec{\hat{r}}) \vec{\hat{r}} -\vec{m} \right) +
%    \frac{2\mu_0}{3} \vec{m} \delta^3(\vec{r})

switch length(r),

  case 2, % 2D case

	  [X, Y] = deal(r{:});
    X = X - rm(1);
    Y = Y - rm(2);
		R2 = X.^2 + Y.^2;
		R = sqrt(R2);
		mCrossRu2 = (m(1)*Y - m(2)*X).^2./R2;
		A = mCrossRu2./R;

  case 3, % 3D case

    [X, Y, Z] = deal(r{:});
    X = X-rm(1);
    Y = Y-rm(2);
    Z = Z-rm(3);
		R2 = X.^2 + Y.^2 + Z.^2;
		R3 = R2.^(3/2);
		mCrossRu2 = (m(1)*Y - m(2)*X).^2./R2;
		A = mCrossRu2./R;

	otherwise

	  error('length(r) should either be 2 or 3');

end


