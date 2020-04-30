Temporal repetition detection for ground visibility assessment
==============================================================

Version 0.7 - April 30, 2020
by rafael grompone von gioi <grompone@gmail.com>


Files
-----

README.txt     - This file
COPYING        - GNU AFFERO GENERAL PUBLIC LICENSE Version 3
Makefile       - Compilation instructions for 'make'
visibility.c   - ANSI C89 code for the visibility detector
iio.c          - Image input/output routines from https://github.com/mnhrdt/iio
iio.h          - Image input/output routines from https://github.com/mnhrdt/iio
data           - Set of ten test images


Compiling
---------

The compiling instruction is just

  make

from the directory where the source codes and the Makefile are located.


Test
----

To verify a correct compilation you can apply the algorithm to the test images
provided with:

  make test


Running
-------

A typical execution is:

  ./visibility 500 image0.tif image1.tif image2.tif image4.tif image5.tif

where 500 is the selected value for the parameter 'lambda' (the size of the
grain filter regions to remove in pixels) and image0.tif to image5.tif are
the satellite time series.  The output of the method is provided as a set
of PNG images with number from 0 to N as 000.png, 001.png, 002.png and so on.
Each of the output visibility mask contain 0 or 255 pixel values: zero
corresponds to visible regions while 255 to an non-visible one.

The input images can be in any image format handled by Enric Meinhardt's
IIO library (https://github.com/mnhrdt/iio).


Copyright and License
---------------------

Copyright (c) 2019-2020 rafael grompone von gioi <grompone@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
