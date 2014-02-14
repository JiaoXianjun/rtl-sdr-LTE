function d=lte_conv_encode(c);

% d=lte_conv_encode(c);
%
% Tail biting convolutional encoder according to LTE specs.

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
error(chk_param(c,'c','vector','real','integer','>=',0,'<=',1,'length>=',6));

n_bits=length(c);

d=NaN(3,n_bits);
s=c(end:-1:end-5);
for t=1:n_bits
  d(1,t)=xor(c(t),xor(s(2),xor(s(3),xor(s(5),s(6)))));
  d(2,t)=xor(c(t),xor(s(1),xor(s(2),xor(s(3),s(6)))));
  d(3,t)=xor(c(t),xor(s(1),xor(s(2),xor(s(4),s(6)))));
  s=[c(t) s(1:5)];
end
%s

