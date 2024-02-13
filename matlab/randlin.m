function r = randlin(X1,X2,varargin)
% function r = randlin(X1,X2,varargin)

%
% $Id: randlin.m,v 1.2 2011/03/26 15:36:05 patrick Exp $
%
% Copyright (c) 2009-2011  Patrick Guio <patrick.guio@gmail.com>
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

ru = rand(varargin{:});

% y = ax + b 
% where 
% a = (y2-y1)/(x2-x1)
% b = (y1*x2-y2*x1)
a = (X2(2)-X1(2))/(X2(1)-X1(1));
b = (X1(2)*X2(1)-X2(2)*X1(1))/(X2(1)-X1(1));

if 0
S = (X2(1)-X1(1))*(a/2*(X2(1)+X1(1))+b)
x = linspace(X1(1), X2(1), 10);
plot(x, a*x+b)
axis equal
pause
end

if a~=0,

  A = a/2;
  B = b;
  C = -a/2*(X1(1)^2+ru*(X2(1)^2-X1(1)^2)) - b*(X1(1)+ru*(X2(1)-X1(1)));

  delta = B^2-4*A*C;

  if any(delta<0),
    error('delta should be positive')
  end

  if 0
    plot(ru,delta)
  end

  r = (-B+sqrt(delta))./(2*A);

else

 r = (X2(1)-X1(1))*ru+X1(1);

end
