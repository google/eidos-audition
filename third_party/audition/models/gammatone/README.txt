1. Original implementation by Malcolm Slaney:
---------------------------------------------
Malcolm Slaney (1998) "Auditory Toolbox Version 2", Technical Report
#1998-010, Interval Research Corporation, 1998:
https://engineering.purdue.edu/~malcolm/interval/1998-010

Original C++ version by João Felipe Santos:
-------------------------------------------
https://github.com/jfsantos/AuditoryFilterbanks

The MIT License (MIT)

Copyright (c) 2015 João Felipe Santos

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

2. The gammatone version from Ning Ma (University of Sheffield):
----------------------------------------------------------------
http://staffwww.dcs.shef.ac.uk/people/N.Ma/resources/gammatone/

*=========================================================================
 * An efficient C implementation of the 4th order gammatone filter
 *-------------------------------------------------------------------------
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *-------------------------------------------------------------------------
%
%  [bm, env, instp, instf] = gammatone_c(x, fs, cf, hrect)
%
%  x     - input signal
%  fs    - sampling frequency (Hz)
%  cf    - centre frequency of the filter (Hz)
%  hrect - half-wave rectifying if hrect = 1 (default 0)
%
%  bm    - basilar membrane displacement
%  env   - instantaneous envelope
%  instp - instantaneous phase (unwrapped radian)
%  instf - instantaneous frequency (Hz)
%
%
%  The gammatone filter is commonly used in models of the auditory system.
%  The algorithm is based on Martin Cooke's Ph.D work (Cooke, 1993) using
%  the base-band impulse invariant transformation. This implementation is
%  highly efficient in that a mathematical rearrangement is used to
%  significantly reduce the cost of computing complex exponentials. For
%  more detail on this implementation see
%  http://www.dcs.shef.ac.uk/~ning/resources/gammatone/
%
%  Once compiled in Matlab this C function can be used as a standard
%  Matlab function:
%  >> mex gammatone_c.c
%  >> bm = gammatone_c(x, 16000, 200);
%
%  Ning Ma, University of Sheffield
%  n.ma@dcs.shef.ac.uk, 09 Mar 2006
%
 * CHANGES:
 * 2012-05-30 Ning Ma <n.ma@dcs.shef.ac.uk>
 *   Fixed a typo in the implementation (a5). The typo does not make a lot
 *   of difference to the response. Thanks to Vijay Parsa for reporting
 *   the problem.
 *
 * 2010-02-01 Ning Ma <n.ma@dcs.shef.ac.uk>
 *   Clip very small filter coefficients to zero in order to prevent
 *   gradual underflow. Arithmetic operations may become very slow with
 *   subnormal numbers (those smaller than the minimum positive normal
 *   value, 2.225e-308 in double precision). This could happen if the
 *   input signal cotains many zeros (e.g. impulse responses). Thanks to
 *   John Culling for reporting the problem.
