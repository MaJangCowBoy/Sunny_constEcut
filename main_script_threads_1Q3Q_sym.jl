#? execute like julia_wh --threads=104 , etc

#? tested with Sunny v0.7.4

using Sunny, HDF5, ProgressBars, CairoMakie
using LinearAlgebra, Statistics, Rotations, Printf
include("function_bundle.jl");

#? number of threads
npar = Threads.nthreads();
#? number of threads

#? Data loading part ?#
data = Matrix{Float64}(undef, 0, 6);
f = open("LT_minimize_1p3000_trim.dat", "r");
for line in eachline(f)  
  tok = parse.(Float64,split(line));  global data = [data; tok'];
end
close(f);  #* data stored like, j2 - jc1 - jc2 - Q1 - Q2 - Q3.
#? Data loading part ?#

#? Basic parameters ?#
kernel = gaussian(fwhm=0.2);  formfactors = [1 => FormFactor("Co2")];
cryst = Crystal("CoTaS.cif",symprec=1e-3);  CoTa3S6 = subcrystal(cryst, "Co");
J1 = 1.311;  Kz = -0.001;
b1 = [0.00, 0.03, 0.06];  B1 = J1 .* b1;
#? Basic parameters ?#

#? define q-points for LSWT ?#
axis1 = [ 1.0, 0.0, 0.0];  N1 = 240;  axis2 = [-0.5, 1.0, 0.0];  N2 = 240;
qgridA, range1, range2, norm1, norm2 = define_qgrid(cryst,axis1,axis2,N1,N2);
qpath1, qpath2, range1, range2, norm1, norm2 = define_qline(cryst,axis1,axis2,N1,N2);
axis1 = [ 0.0, 1.0, 0.0];  N1 = 240;  axis2 = [-1.0, 0.5, 0.0];  N2 = 240;
qgridB, range1, range2, norm1, norm2 = define_qgrid(cryst,axis1,axis2,N1,N2);
axis1 = [-1.0, 1.0, 0.0];  N1 = 240;  axis2 = [-0.5,-0.5, 0.0];  N2 = 240;
qgridC, range1, range2, norm1, norm2 = define_qgrid(cryst,axis1,axis2,N1,N2);
#? define q-points for LSWT ?#