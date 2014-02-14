function x=chi2inv(p,k)

% x=chi2inv(p,k)
%
% Returns the value x such that chi2cdf(x,k) is p.
%
% This function is very slow, but works...

% Copyright 2012 Evrytania LLC (http://www.evrytania.com)
%
% Written by James Peroulas <james@evrytania.com>
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

error(nargchk(2,2,nargin));
error(chk_param(p,'p','scalar','real','>=',0,'<',1));
error(chk_param(k,'k','scalar','real','integer','>',0));

try_lo=0;
try_hi=100;
% First, search for an appropriate try_hi
while (chi2cdf(try_hi,k)<p)
  try_hi=try_hi^2;
end

% Now, search for the solution
while (try_hi-try_lo>10*eps(try_hi))
  try_mid=mean([try_hi try_lo]);
  if (chi2cdf(try_mid,k)>p)
    try_hi=try_mid;
  else
    try_lo=try_mid;
  end
end

x=try_hi;

