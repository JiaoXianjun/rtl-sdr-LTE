function f=chi2pdf(x,k)

% f=chi2pdf(x,k)
%
% Return the pdf of a chi-squared distribution with k degrees of freedom
% evaluated at x.

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
error(chk_param(x,'x','numeric','real'));
error(chk_param(k,'k','scalar','real','integer','>=',1));

f=NaN(size(x));
f(x<0)=0;
xgt0=x>=0;
f(xgt0)=x(xgt0).^(k/2-1).*exp(-x(xgt0)/2)/(2^(k/2)*gamma(k/2));

