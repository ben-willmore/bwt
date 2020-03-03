function [si, cx] = bwt_v1(im)
% function decomp = bwt(im)
%
% BWT decomposition into idealised V1-like units.
% Version 1.2
%
% Arguments:
%  im: A square image; side length must be a power of 3
%
% Result:
%  si: The half-wave rectified simple cell decomposition.
%  si is a cell array of cell arrays such that:
%  si{1} is the finest-resolution decomposition (n/3 x n/3), ... si{end} is the coarsest (1x1)
%  si{n} contains cell arrays containing the decomposition for a single orientation and phase.
%  si{n} is ordered: si{n}{1}: or=0;ph=0
%                    si{n}{2}: or=0;ph=90
%                    si{n}{3}: or=0;ph=180
%                    si{n}{4}: or=0;ph=270
%                    si{n}{5}: or=45;ph=0, ... up to or=135;ph=270
%  Each si{n} has 16 elements in total
%
%  cx: The energy model complex cell decomposition
%  cx is a cell array of cell arrays such that:
%  cx{1} is the finest-resolution decomposition (n/3 x n/3), ... cx{end} is the coarsest (1x1)
%  cx{n} is ordered: cx{n}{1}: or=0
%                    cx{n}{2}: or=45
%                    cx{n}{3}: or=90
%                    cx{n}{4}: or=135
%  Each cx{n} has 4 elements in total
%
% Citation:
%  Willmore B, Prenger RJ, Wu MC and Gallant JL (2008). The Berkeley
%  Wavelet Transform: A biologically-inspired orthogonal wavelet transform.
%  Neural Computation 20:6, 1537-1564
%
% The article is available at:
%  <http://dx.doi.org/10.1162/neco.2007.05-07-513>
%
% Copyright (c) 2020 Ben Willmore
%
% Permission is hereby granted, free of charge, to any person
% obtaining a copy of this software and associated documentation
% files (the "Software"), to deal in the Software without
% restriction, including without limitation the rights to use,
% copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the
% Software is furnished to do so, subject to the following
% conditions:
%
% The above copyright notice and this permission notice shall be
% included in all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
% HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
% WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
% FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
% OTHER DEALINGS IN THE SOFTWARE.

sz = size(im);

if (length(sz) ~= 2) || (sz(1) ~= sz(2))
  disp('Input must be square');
  decomp = nan;
  return;
end

sz = sz(1);

numlevels = log(sz)/log(3);

if ( (numlevels-floor(numlevels)) > abs(numlevels)*eps )
  disp('Input side length must be a power of 3');
  decomp = nan;
  return;
end

decomp = zeros(sz);
si = {};
cx = {};

for level = 1:numlevels
  decomp_thislevel = bwt_onelevel(im);
  ssz = size(decomp_thislevel,1);
  decomp(end-ssz+1:end,1:ssz) = decomp_thislevel/((3^level)^2);

  % construct half wave rectified simple cells
  %  [ 90-odd  135-odd  45-even ]
  %  [ 90-even 135-even 45-odd  ]
  %  [ DC      0-even   0 -odd  ]
  sections = {};
  for yy = 1:3
    for xx = 1:3
      sections{end+1} = decomp_thislevel((yy-1)*ssz/3+1:yy*ssz/3, (xx-1)*ssz/3+1:xx*ssz/3);
    end
  end
  si_norect = {sections{9} sections{8} sections{6} sections{3} sections{1} sections{4} sections{2} sections{5}};
  si_rect = {};
  for ii = 1:2:length(si_norect)
    si_rect{end+1} =  max(si_norect{ii}, 0);
    si_rect{end+1} =  max(si_norect{ii+1}, 0);
    si_rect{end+1} = -min(si_norect{ii}, 0);
    si_rect{end+1} = -min(si_norect{ii+1}, 0);
  end
  si{level} = si_rect;
  cx{level} = {si_norect{1}.^2 + si_norect{2}.^2 ...
               si_norect{3}.^2 + si_norect{4}.^2 ...
               si_norect{5}.^2 + si_norect{6}.^2 ...
               si_norect{7}.^2 + si_norect{8}.^2};

  if (ssz>1)
    ssz = ssz/3;
    im = decomp_thislevel(end-ssz+1:end,1:ssz);
  end

end
