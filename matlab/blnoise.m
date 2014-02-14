function y=blnoise(n_samp)

% y=blnoise(n_samp)
%
% Return a signal composed of i.i.d. complex Gaussian samples of length
% n_samp.

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

error(nargchk(1,1,nargin));
error(chk_param(n_samp,'n_samp','scalar','real','integer','>=',0));

y=complex(randn(1,n_samp),randn(1,n_samp))/sqrt(2);

