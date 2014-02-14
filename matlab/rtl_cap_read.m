function c=rtl_cap_read(filename)

% c=rtl_cap_read(filename)
%
% Read a capture file created by rtl_sdr and place the result in to the
% vector c.

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
error(chk_param(filename,'filename','string'));

f=fopen(filename,'rb');
c_raw=fread(f,[2 Inf],'uint8');

c=complex(c_raw(1,:)-127,c_raw(2,:)-127)/128;

