The code in this directory replaces the usual source code that is distributed with the
BOOZ_XFORM sources supplied with STELLOPT. There are five new source files (boozer_metric.f,
write_polcut.f, root.f, vcoords_gijb.f and boozer_gij.f) and small modifications to the standard
BOOZ_XFORM distribution files (boozer_xform.f, booz_params.f, boozer.f, allocate_boozer.f).
These versions have been modified so that when xbooz_xform is run with an
additional command line argument "far" it produces the woutb file used by FAR3d.
woutb is a binary file that contains data about the radial grid, the Fourier modes
used to represent the Boozer transform data, and Fourier amplitudes for each
surface of B, R, Z, Phi, and various metric elements in Boozer coords. as used by FAR3d.
This modification was provided by Luis Garcia of the Universidad Carlos III-Leganes, Spain on 4-3-2018.

To run in the normal way, one types: xbooz_xform input.boz
To run and generate the woutb file: xbooz_xform input.boz far

In addition to the to the woutb file, you get some text diagnostic files with the Boozer
components of R, Z, phi, B,the Jacobian and the metric elements. Also files are provided
with the reconstruction of some flux surfaces from the vmec and boozer coordinates (R,Z,etc.).
