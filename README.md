vector_compression
==================

The `vector_compression` library contains some functions for compressing
and decompressing 3D and 4D (unit) vectors. It was developed for the
participation of Team [NimbRo Rescue] in the [DARPA Robotics Challenge].
The challenge regulations enforce a bandwidth limit of 9600 bit/s on the
communication link with the robot, allowing only small one-second bursts
of high-bandwidth communication. Thus, we had to compress commands and feedback
that we wanted to sent at higher rates (up to 10Hz) as much as possible.

While generic compression (e.g. LZMA) is useful, you can achieve higher
compression ratios and deterministic behavior using compression routines
tailored to the specific task. Such algorithms were often originally developed
for fast data transfer between CPU and GPU, but work equally well in this
low-bandwidth situation.

This small library contains routines for

 * (de)compression of signed floating point numbers to/from arbitrary bit widths
 * (de)compression of Quaternions to/from 5 bytes using the [hemi-oct encoding]
 * (de)compression of 3D vectors to arbitrary bit widths using a face-centered
   cubic packing (FCC) lattice.

It also contains unit tests for most routines.

[Nimbro Rescue]: http://www.ais.uni-bonn.de/nimbro/Rescue/index.html
[DARPA Robotics Challenge]: http://theroboticschallenge.org/
[hemi-oct encoding]: http://jcgt.org/published/0003/02/01/

Building & Usage
----------------

The library can be built as a ROS catkin package (simply include it in your
workspace) or as a plain CMake project:

```
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo ..
make
sudo make install
```

License
-------

`vector_compression` is released under the BSD 3-clause license. It contains
a header file from the ROS `angles` package, which is also released under BSD.

Authors & Contact
-----------------

Author: `Max Schwarz <max.schwarz@uni-bonn.de>` <br />
Please contact me if you have questions!
