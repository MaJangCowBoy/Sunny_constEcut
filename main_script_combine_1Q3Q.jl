#? tested with Sunny v0.7.4

using Sunny, HDF5, ProgressBars, CairoMakie
using LinearAlgebra, Statistics, Rotations, Printf
include("function_bundle.jl");

#? mode selection part ?#
sweep_mode = ARGS[1];  # energies = [1.5];
#? mode selection part ?#

#? Basic parameters ?#
kernel = lorentzian(fwhm=2.0);  formfactors = [1 => FormFactor("Co2")];

J1 = 1.60;  Kz = -0.001;
b1 = [0.00, 0.02, 0.04, 0.06];  B1 = J1 .* b1 ;

j2 = parse(Float64,ARGS[2]);  jc1 = parse(Float64,ARGS[3]);  jc2 = parse(Float64,ARGS[4]);
F = jc1 + jc2;  G = jc1 - jc2 * 0.5;
j3 = 1/2 * (1 - (F*G - jc2*F/2 + 2*jc2*G )/√(F*F-2*F*G+4*G*G));
J2 = j2 * J1;  J3 = j3 * J1;  Jc1 = jc1 * J1;  Jc2 = jc2 * J1;
Jc1mat = Jc1 * ([1 0 0; 0 1 0; 0 0 1] + 0.001 * dmvec([0, 0, 1]));

energies = [0.5 * (J1+J2)]; # meV
  # added for trend check purpose

#? Basic parameters ?#

#? define system part ?#
cryst = Crystal("CoTaS.cif",symprec=1e-3);  CoTa3S6 = subcrystal(cryst, "Co");

if length(B1) == 1
  B1 = B1[1];
  sys = define_system(CoTa3S6,"3Q",J1,B1,J2,J3,Jc1mat,Jc2,Kz,(3,3,1));
  keyword = "3Q";  sys = system_initialize(sys, keyword, J1);
  measure = ssf_perp(sys; formfactors);  swt = SpinWaveTheory(sys; measure);
else
  swt = [];
  for Bi in B1
    sys_i = define_system(CoTa3S6,"3Q",J1,Bi,J2,J3,Jc1mat,Jc2,Kz,(3,3,1));
    keyword = "3Q";  sys_i = system_initialize(sys_i, keyword, J1);
    measure = ssf_perp(sys_i; formfactors);  swt_i = SpinWaveTheory(sys_i; measure);  
    global swt = [swt; swt_i];
  end
end

sys1 = define_system(CoTa3S6,"1Q",J1,B1,J2,J3,Jc1mat,Jc2,Kz,(3,3,1));
keyword = "1Q_1";  sys1 = system_initialize(sys1, keyword, J1);
measure = ssf_perp(sys1; formfactors);  swt1 = SpinWaveTheory(sys1; measure);
sys2 = define_system(CoTa3S6,"1Q",J1,B1,J2,J3,Jc1mat,Jc2,Kz,(3,3,1));
keyword = "1Q_2";  sys2 = system_initialize(sys2, keyword, J1);
measure = ssf_perp(sys2; formfactors);  swt2 = SpinWaveTheory(sys2; measure);
sys3 = define_system(CoTa3S6,"1Q",J1,B1,J2,J3,Jc1mat,Jc2,Kz,(3,3,1));
keyword = "1Q_3";  sys3 = system_initialize(sys3, keyword, J1);
measure = ssf_perp(sys3; formfactors);  swt3 = SpinWaveTheory(sys3; measure);
#? define system part ?#

#? define q-points and calculate intensities ?#
axis1 = [ 1.0, 0.0, 0.0];  N1 = 60;  axis2 = [-0.5, 1.0, 0.0];  N2 = 60;
# if sweep_mode == "2D"
  qgrid, range1, range2, norm1, norm2 = define_qgrid(cryst,axis1,axis2,N1,N2);

  data_3Q = zeros(Float64,length(B1),N1,N2);
  if length(swt) == 1
    res = intensities(swt, qgrid, energies, kernel);
    global data_3Q[1,:,:] = res.data[1,:,:];
  else
    for (i,swt_i) in enumerate(swt)
      res_i = intensities(swt_i, qgrid; energies, kernel);
      global data_3Q[i,:,:] = res_i.data[1,:,:];
    end
  end

  res1 = intensities(swt1, qgrid; energies, kernel);
  res2 = intensities(swt2, qgrid; energies, kernel);
  res3 = intensities(swt3, qgrid; energies, kernel);
  data_1Q = res1.data[1,:,:] + res2.data[1,:,:] + res3.data[1,:,:];

# elseif sweep_mode == "1D"
  qpath1, qpath2, range1, range2, norm1, norm2 = define_qline(cryst,axis1,axis2,N1,N2);

  data_1_3Q = zeros(Float64,length(B1),N1);  data_2_3Q = zeros(Float64,length(B1),N2);
  if length(swt) == 1
    res1 = intensities(swt, qpath1; energies, kernel);
    global data_1_3Q[1,:] = res1.data[1,:];
    res2 = intensities(swt, qpath2; energies, kernel);
    global data_2_3Q[1,:] = res2.data[1,:];
  else
    for (i,swt_i) in enumerate(swt)
      res1_i = intensities(swt_i, qpath1; energies, kernel);
      global data_1_3Q[i,:] = res1_i.data[1,:];
      res2_i = intensities(swt_i, qpath2; energies, kernel);
      global data_2_3Q[i,:] = res2_i.data[1,:];
    end
  end

  res1_1 = intensities(swt1, qpath1; energies, kernel);
  res2_1 = intensities(swt2, qpath1; energies, kernel);
  res3_1 = intensities(swt3, qpath1; energies, kernel);
  data_1_1Q = res1_1.data[1,:] + res2_1.data[1,:] + res3_1.data[1,:];

  res1_2 = intensities(swt1, qpath2; energies, kernel);
  res2_2 = intensities(swt2, qpath2; energies, kernel);
  res3_2 = intensities(swt3, qpath2; energies, kernel);
  data_2_1Q = res1_2.data[1,:] + res2_2.data[1,:] + res3_2.data[1,:];
# end
#? define q-points and calculate intensities ?#


#? data saving part ?#
if length(B1) == 1
  tail    = @sprintf("B1_%+.2f_J2_%+.2f_Jc1_%+.2f_Jc2_%+.2f",B1,j2,jc1,jc2);
  tail    = replace(tail,"." => "p","-" => "M","+" => "P");
else
  tail    = @sprintf("B1_%+.2f_%+.2f_J2_%+.2f_Jc1_%+.2f_Jc2_%+.2f",B1[1],B1[end],j2,jc1,jc2);
  tail    = replace(tail,"." => "p","-" => "M","+" => "P");
end
figname = @sprintf("data_figure_combine_1Q3Q_%s_%s.png",sweep_mode,tail);
h5name  = @sprintf("data_h5file_combine_1Q3Q_%s_%s.h5", sweep_mode,tail);

if sweep_mode == "2D"
  # fig = Figure();  ax = Axis(fig[1, 1], aspect = norm1/norm2);
  # heatmap!(ax, range1, range2, data[:,:];  colormap = (:viridis, 1), colorrange = (1,100));
  # fig, ax = add_BZ_boundary(fig, ax);  xlims!(ax, -1, 1);  ylims!(ax, -√3/2, √3/2);  
  # save(figname,fig);
  export_to_h5file2D_combine_1Q3Q(h5name,data_1Q,data_3Q,range1,range2,norm1,norm2,b1,j2,jc1,jc2);
elseif sweep_mode == "1D"
  # fig = Figure();  ax1 = Axis(fig[1, 1]);  ax2 = Axis(fig[2, 1]);
  # plot!(ax1, range1, data_1; color = :blue, linewidth = 2);
  # plot!(ax2, range2, data_2; color = :red, linewidth = 2);
  # save(figname,fig);
  export_to_h5file1D_combine_1Q3Q(h5name,data_1_1Q,data_2_1Q,data_1_3Q,data_2_3Q,range1,range2,norm1,norm2,b1,j2,jc1,jc2);
end

export_to_h5file_combine_1Q3Q
#? data saving part ?#