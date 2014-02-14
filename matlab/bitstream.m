function b=bitstream(l,p)

% b=bitstream(l,p)
%
% Return a stream of l bits where the probability of a '1' is p.
% p defaults to 0.5.

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

error(nargchk(1,2,nargin));
error(chk_param(l,'l','scalar','integer','real','>=',0));
if (nargin<2)
  p=0.5;
end
error(chk_param(p,'p','scalar','real','>=',0,'<=',1));

b=double(rand(1,l)<p);

