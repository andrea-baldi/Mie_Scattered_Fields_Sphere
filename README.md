# Mie_Scattered_Fields_Sphere
Plot the fields scattered by a spherical particle using Mie theory

This script reads the relative permittivity data of a material and uses Mie theory to compute the scattering efficiency for a spherical particle of that material embedded in a lossless medium. The permittivity file has to be a tab-delimited text file with three columns: energy (in eV), epsilon1, epsilon2.

The code then lets you manually choose at which energy you want to map the scattered fields, which are calculated using the equations of chapter 4 in Bohren and Huffman. In agreement with the book's system of reference, the incident plane wave propagates along the z direction, with the electric
field polarized along the x-axis. The polar angle 'theta', defined as the angle between the scattered vector and the z-axis, spans from 0 degrees (forward scattering) to 180 degrees (back scattering), the azimuthal angle 'phi', defined as the angle between the x-axis and the projection of the scattered vector on the xy-plane, spans from 0 degrees to 360 degrees. The code needs the function "pin_andrea.m", corresponding to eq. 4.47.

Finally, you have the option to plot the near-field spectrum at a specific location of the map. In most cases the near-field spectrum at a given location (Figure 3) is very different from the scattering efficiency of the overall particle (Figure 1).
