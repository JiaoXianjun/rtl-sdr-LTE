function t=Tc(v,fc)

% t=Tc(v,fc)
%
% Return the channel coherence time for a mobile moving v m/s and communicating
% over carrier frequenc fc in Hz.

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
error(chk_param(v,'v','scalar','real','>=',0));
error(chk_param(fc,'fc','scalar','real','>=',0));

dop_spread=2*doppler_shift(fc,v);
t=1/4/dop_spread;

