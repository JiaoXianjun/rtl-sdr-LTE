function syms=qam(bits,mod,map)

% syms=qam(bits,mod,map)
%
% Convert 'bits' into symbols using modulation 'mod'.
% Supported modulations are 'QAM','QAM16','QAM64','QAM256'
%
% 'map' indicates how bits are mapped to symbols. For QAM, 'map'
% may be [0 1 ; 2 3] which indicates that the bitsequence 00 maps
% to the upper left symbol and 10 maps to the lower left symbol.
%
% 'map' can also be 'LTE' (default) which produces the LTE mapping
% between bits and symbols.
%
% Output power is normalized to 0dB.

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

error(nargchk(2,3,nargin));
error(chk_param(bits,'bits','vector','real','integer','>=',0,'<=',1));
error(chk_param(mod,'string'));
if (nargin<3)
  map='LTE';
end
if (~ischar(map))
  error(chk_param(map,'map','real','integer','>=',0));
end

% Second order checking and derivations
n_bits=length(bits);
switch (lower(mod))
  case 'qam'
    bps=2;
  case 'qam16'
    bps=4;
  case 'qam64'
    bps=6;
  case 'qam256'
    bps=8;
  otherwise
    error('Unrecognized modulation specified.');
end
map_dim=2^(bps/2);
if (ischar(map))
  if (strcmpi(map,'LTE'))
    switch (bps)
      case 2
        map=[2 0 ; 3 1];
      case 4
        map=[11 9 1 3; 10 8 0 2; 14 12 4 6; 15 13 5 7];
      case 6
        map=[47 45 37 39 7 5 13 15; 46 44 36 38 6 4 12 14; 42 40 32 34 2 0 8 10; 43 41 33 35 3 1 9 11; 59 57 49 51 19 17 25 27; 58 56 48 50 18 16 24 26; 62 60 52 54 22 20 28 30; 63 61 53 55 23 21 29 31];
      otherwise
        error(['no lte mapping is known for ' mod]);
    end
  else
    error('unrecognized mapping specified.');
  end
else
  error(chk_param(map,'map','integer','real','>=',0,'<=',2^bps-1));
end
if ((length(map(:))~=2^bps)||(size(map,1)~=map_dim)||size(map,2)~=map_dim)
  error(sprintf('map matrix dimensions must be %ix%i',map_dim,map_dim));
end

% Final checking
if (length(bits)/bps~=floor(length(bits)/bps))
  error(sprintf('bitstream length must be a multiple of %i',bps));
end
if (length(unique(map(:)))~=2^bps)
  error('map matrix contains duplicate entries');
end

% Create the mapping between integers and complex modulation symbols.
const=complex(repmat(1:map_dim,map_dim,1),repmat(transpose(map_dim:-1:1),1,map_dim));
const=const-mean(const(:));
const=const/sqrt(sigpower(const(:)));
map_flat=NaN(1,2^bps);
for t=1:map_dim
  for m=1:map_dim
    map_flat(map(t,m)+1)=const(t,m);
  end
end
%sprintf('%i ',round(imag((map_flat)*sqrt(42))))

% Combine bits into integers
int=zeros(1,n_bits/bps);
for t=1:bps
  int=int*2+bits(t:bps:end);
end

% Create the symbols
syms=map_flat(int+1);

