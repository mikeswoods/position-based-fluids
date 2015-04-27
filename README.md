# Real(ish)-time position based fluid simulation

##### CIS 563: Physically Based Animation, Spring 2015
##### By [Michael Woods](mailto:micwoods@seas.upenn.edu) & [Michael O'Meara](mailto:momeara@seas.upenn.edu)
###### Instructor: Ladislav Kavan
###### TAs: Yu Wang, Harmony Li, Xinjie Ma, Ying Li, and Tiantian Liu

Background
----------

For the purpose of this project, we had two main goals in mind: to correctly and _cleanly_
implement a real-time, particle-based fluid simulation as outlined in the 2013 paper 
[Position Based Fluids](http://mmacklin.com/pbf_sig_preprint.pdf) by Matthias Müller and Miles Macklin.
To achieve the goal of real-time simulation, we exploited the parallel processing 
power of the GPU as much as possible, either through the use of OpenGL for
direct rendering, or OpenCL for generalized parallel computation.

Implementation details
----------------------
For the most part, every major feature outlined by Müller and Macklin was 
implemented in our simulator:
- Particle incompressibility
- Simulated surface tension using an artificial pressure term
- Vorticity confinement to mitigate numerical dissipation issues
- XSPH viscosity

Enhancements 
------------
- All major components of the simulation were implemented as OpenCL kernels,
  allowing the code to run directly on the GPU
- As a potential improvement over the fast nearest-neighbor finding method 
  employed by Müller and Macklin [(Green 2008)](http://developer.download.nvidia.com/presentations/2008/GDC/GDC08_ParticleFluids.pdf), we implemented the parallel counting sort technique as described in [FAST FIXED-RADIUS NEAREST NEIGHBORS: INTERACTIVE MILLION-PARTICLE FLUIDS](http://on-demand.gputechconf.com/gtc/2014/presentations/S4117-fast-fixed-radius-nearest-neighbor-gpu.pdf) by Hoetzlein/NVIDIA (2013)
- The bounds of the simulation can be dynamically changed at runtime according
  to a number of predefined animation patterns, such as sine wave, sawtooth wave,
  and one-shot compression. This allows us to create a number of novel fluid
  motion effects like crashing waves, etc. with ease on-the-fly

Tooling, third Party libraries, and frameworks
----------------------------------------------

- [openFrameworks](http://www.openframeworks.cc/), community [contributed](http://www.openframeworks.cc/list-info)
  * A great, generalized C++ framework for creative coding and experimentation
- [ofxMSAOpenCL](https://github.com/memo/ofxMSAOpenCL), Memo Akten
  * A openFrameworks addon wrapper that removes much of the boilerplate code 
    associated with OpenCL

References
----------

- [Position Based Fluids](http://mmacklin.com/pbf_sig_preprint.pdf)
  * Matthias Müller and Miles Macklin (2013)
- [Position Based Fluids, SIGGRAPH 2013 slides](http://mmacklin.com/pbf_slides.pdf)
  * Matthias Müller and Miles Macklin
- [Particle-based Fluid Simulation based Fluid Simulation](http://developer.download.nvidia.com/presentations/2008/GDC/GDC08_ParticleFluids.pdf)
  * Simon Green (2008)
- [FAST FIXED-RADIUS NEAREST NEIGHBORS: INTERACTIVE MILLION-PARTICLE FLUIDS](http://on-demand.gputechconf.com/gtc/2014/presentations/S4117-fast-fixed-radius-nearest-neighbor-gpu.pdf)
  * Rama C. Hoetzlein, Graphics Devtech, NVIDIA (2013)
- [OpenCL Parallel Prefix Sum (aka Scan) Example](https://developer.apple.com/library/mac/samplecode/OpenCL_Parallel_Prefix_Sum_Example/Introduction/Intro.html)
  * Apple Inc. Updated: 2009-09-23
- [Point sprites as spheres in OpenGL3.3](http://mmmovania.blogspot.de/2011/01/point-sprites-as-spheres-in-opengl33.html)
  * Movania Muhammad Mobeen (2011)
- [Scaling point sprites with distance](http://gamedev.stackexchange.com/a/54492)
- [Allocating Local Memory](http://www.openclblog.com/2014/10/allocating-local-memory.html)
  * Matt Scarpino (2014)
- [OpenCL Sorting](http://www.bealto.com/gpu-sorting.html)
  * Eric Bainville, June 2011

Building 
--------

#### Hardware prerequisites
A graphics card that supports OpenCL v1.1+ and OpenGL 3+ (a.k.a. "modern OpenGL")

#### Software prerequisites

##### OSX

- XCode 6.x

- OSX 10.9 and higher, although it might compile on older versions. Unfortunately,
we are unable to test for compatibility on older versions of OSX.

##### Windows 7+

- Visual Studio 2012 Professional (although the Basic edition may suffice)

- The necessary OpenCL runtime/SDK provided by the manufacturer of the
  the graphics card/devices:

  - For __Intel__ cards and chipsets:
    1. Download the latest Intel graphics [drivers](https://software.intel.com/en-us/articles/opencl-drivers)
    2. Download the free edition of [OpenCL™ Code Builder with Intel® INDE](https://software.intel.com/en-us/intel-opencl) 
       (formerly Intel® SDK for OpenCL™ Applications). This is needed for `OpenCL.lib` for Visual Studio to link against.
       If it is installed correctly, you should see the environment variable __`INTELOCLSDKROOT`__ defined under
       `System > Advanced System Settings > Environment Variable > System variables`

  - For __Nvidia__ cards and chipsets: 
    1. *TODO*

  - For __AMD__ cards and chipsets:
    1. *TODO*

  - For __ATI__ cards and chipsets: 
    1. *TODO*
