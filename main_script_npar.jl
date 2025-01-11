#? Please execute with "julia_wh --threads=auto"
#? or "julia_wh --threads=4" etc.
#? tested with Sunny v0.7.4

using Sunny, HDF5, ProgressBars, CairoMakie
using LinearAlgebra, Statistics, Rotations
include("function_bundle.jl");

#? mode selection part
main_keyword = "1Q";  sweep_mode = "2D";
#? mode selection part

npar = Threads.nthreads();

filename = "LT_minimize_0p0000.dat";
Js_arr = Matrix{Float64}(undef,0,3);

file = open(filename,"r");
for row in eachline(file)
  tmp = parse.(Float64,split(row));
  j2 = tmp[1];  jc1 = tmp[2];  jc2 = tmp[3];
  q1 = tmp[4];  q2  = tmp[5];  q3  = tmp[6];

  if q1 > 0.0 && q2 == 0 # if ordering is q00
    Js_arr = [Js_arr; j2 jc1 jc2];
  else    end
end;  close(file);

LenJs = size(Js_arr,1);

Threads.@threads for id in 1:npar
  SttPt = 1 + (id-1) * div(LenJs,npar);
  EndPt = id * div(LenJs,npar);
  if id == npar  EndPt = LenJs;  end

  J1 = 1.6;  B1 = 0.00 * J1;  Kz = -0.001;
  for it1 = SttPt:EndPt
    j2 = Js_arr[it1,1];  jc1 = Js_arr[it1,2];  jc2 = Js_arr[it1,3];
    F = jc1 + jc2;  G = jc1 - jc2 * 0.5;
    j3 = 1/2 * (1 - (F*G - jc2*F/2 + 2*jc2*G )/√(F*F-2*F*G+4*G*G));
    J2 = j2 * J1;  J3 = j3 * J1;  Jc1 = jc1 * J1;  Jc2 = jc2 * J1;
    Jc1mat = Jc1 * ([1 0 0; 0 1 0; 0 0 0] + 0.001 * dmvec([0, 0, 1]));
    
    cryst = Crystal("CoTaS.cif",symprec=1e-3);
    CoTa3S6 = subcrystal(cryst, "Co1");
    sys = System(CoTa3S6, [1 => Moment(s=3/2, g=2)], :dipole)

    set_pair_coupling!(sys, (Si, Sj) -> Si'*J1*Sj + B1*(Si'*Sj)^2, Bond(1, 1, [1, 0, 0]));
    set_exchange!(sys, J2, Bond(1,1,[1,-1,0]));
    set_exchange!(sys, J3, Bond(1,1,[2,0,0]));
    set_exchange!(sys, Jc1mat, Bond(1,2,[0,0,0]));
    set_exchange!(sys, Jc2, Bond(1,2,[1,1,0]));
    set_onsite_coupling!(sys, S -> Kz*S[3]^2, 1);

    sys = repeat_periodically(sys, (3,3,1));
    formfactors = [1 => FormFactor("Co2")];
    kernel = lorentzian(fwhm=2.0);

    axis1 = [ 1.0, 0.0, 0.0];  N1 = 60;  axis2 = [-0.5, 1.0, 0.0];  N2 = 60;
    if sweep_mode == "2D"
      qgrids, range1, range2, norm1, norm2 = define_qgrid(cryst,axis1,axis2,N1,N2);
    elseif sweep_mode == "1D"
      qpath1, qpath2, range1, range2, norm1, norm2 = define_qline(cryst,axis1,axis2,N1,N2);
    end

    if main_keyword == "3Q"
      keyword = "3Q";  sys = system_initialize(sys, keyword);
      dt = 0.1/(J1 * 1.5 * 1.5);  damping = 0.1;  
      langevin = Langevin(dt; damping, kT = 0.0);
      langevin.kT = 0.1 * meV_per_K;  
      for _ in 1:100000
        step!(sys, langevin)
      end
      langevin.kT = 0.0;  
      for _ in 1:100000
        step!(sys, langevin)
      end
      # print_wrapped_intensities(sys);  # plot_spins(sys);
      measure = ssf_perp(sys; formfactors);  swt = SpinWaveTheory(sys; measure);

      if sweep_mode == "2D"
        for qgrid in ProgressBar(qgrids)
          res = intensities(swt, qgrid; energies = [1.5], kernel, kT= 0*meV_per_K);  push!(ress, res);
        end
        data = zeros(size(res.data,2),size(res.data,3))
        for i in axes(res.data, 2), j in axes(res.data, 3)
          data[i,j] = res.data[1,i,j]
        end
      elseif sweep_mode == "1D"
        res1 = intensities(swt, qpath1; energies = [1.5], kernel, kT= 0*meV_per_K);
        res2 = intensities(swt, qpath2; energies = [1.5], kernel, kT= 0*meV_per_K);
        data1 = Vector{Float64}(undef,0);  for i in axes(res1.data, 2)  push!(data1, res1.data[1,i]);  end
        data2 = Vector{Float64}(undef,0);  for i in axes(res2.data, 2)  push!(data2, res2.data[1,i]);  end
        data = [data1, data2];
      end

    elseif main_keyword == "1Q"
      sys1 = sys;  sys2 = sys;  sys3 = sys;
      keyword = "1Q_1";  sys1 = system_initialize(sys3, keyword);
      keyword = "1Q_2";  sys2 = system_initialize(sys3, keyword);
      keyword = "1Q_3";  sys3 = system_initialize(sys3, keyword);
      dt = 0.1/(J1 * 1.5 * 1.5);  damping = 0.1;  
      langevin = Langevin(dt; damping, kT = 0.0);
      langevin.kT = 0.1 * meV_per_K;  
      for _ in 1:100000  
        step!(sys1, langevin);  step!(sys2, langevin);  step!(sys3, langevin);
      end
      langevin.kT = 0.0;  
      for _ in 1:100000  
        step!(sys1, langevin);  step!(sys2, langevin);  step!(sys3, langevin);
      end
      # print_wrapped_intensities(sys);  # plot_spins(sys);
      measure = ssf_perp(sys1; formfactors);  swt1 = SpinWaveTheory(sys1; measure);
      measure = ssf_perp(sys2; formfactors);  swt2 = SpinWaveTheory(sys2; measure);
      measure = ssf_perp(sys3; formfactors);  swt3 = SpinWaveTheory(sys3; measure);

      if sweep_mode == "2D"
        for qgrid in ProgressBar(qgrids)
          res1 = intensities(swt1, qgrid; energies = [1.5], kernel, kT= 0*meV_per_K);  push!(ress1, res1);
          res2 = intensities(swt2, qgrid; energies = [1.5], kernel, kT= 0*meV_per_K);  push!(ress2, res2);
          res3 = intensities(swt3, qgrid; energies = [1.5], kernel, kT= 0*meV_per_K);  push!(ress3, res3);
        end
        res = res1[1];
        for a in axes(res.data, 1), b in axes(res.data, 2), c in axes(res.data, 3)
          res.data[a, b, c] = 0.0;
        end;  for i in axes(ress1, 1)  res.data[:] += res1[i].data[:] + res2[i].data[:] + res3[i].data[:];  end

        data = zeros(size(res.data,2),size(res.data,3))
        for i in axes(res.data, 2), j in axes(res.data, 3)
          data[i,j] = res.data[1,i,j]
        end
      elseif sweep_mode == "1D"
        for qgrid in ProgressBar(qgrids)
          res11 = intensities(swt1, qpath1; energies = [1.5], kernel, kT= 0*meV_per_K);  push!(ress1, res11);
          res12 = intensities(swt1, qpath2; energies = [1.5], kernel, kT= 0*meV_per_K);  push!(ress1, res12);
          res21 = intensities(swt2, qgrid1; energies = [1.5], kernel, kT= 0*meV_per_K);  push!(ress2, res21);
          res22 = intensities(swt2, qgrid2; energies = [1.5], kernel, kT= 0*meV_per_K);  push!(ress2, res22);
          res31 = intensities(swt3, qgrid1; energies = [1.5], kernel, kT= 0*meV_per_K);  push!(ress3, res31);
          res32 = intensities(swt3, qgrid2; energies = [1.5], kernel, kT= 0*meV_per_K);  push!(ress3, res32);
        end
        res1 = res11[1];  res2 = res12[1];
        for a in axes(res1.data, 1), b in axes(res1.data, 2), c in axes(res1.data, 3)
          res1.data[a, b, c] = 0.0;
        end;  for i in axes(ress1, 1)  res1.data[:] += res11[i].data[:] + res21[i].data[:] + res31[i].data[:];  end
        for a in axes(res2.data, 1), b in axes(res2.data, 2), c in axes(res2.data, 3)
          res2.data[a, b, c] = 0.0;
        end;  for i in axes(ress1, 1)  res2.data[:] += res12[i].data[:] + res22[i].data[:] + res32[i].data[:];  end
      
        data1 = Vector{Float64}(undef,0);  for i in axes(res1.data, 2)  push!(data1, res1.data[1,i]);  end
        data2 = Vector{Float64}(undef,0);  for i in axes(res2.data, 2)  push!(data2, res2.data[1,i]);  end
      end
    end

    figname = "figure_"*replace("$(j2)_$(jc1)_$(jc2)","." => "p")*".png";

    if sweep_mode == "2D"
      fig = Figure();  ax = Axis(fig[1, 1], aspect = norm1/norm2);
      heatmap!(ax, range1, range2, data[:,:];  colormap = (:viridis, 1), colorrange = (1,100));
      fig, ax = add_BZ_boundary(fig, ax);
      xlims!(ax, -1, 1);  ylims!(ax, -√3/2, √3/2);  save(figname,fig);
    elseif sweep_mode == "1D"
      fig = Figure();  ax1 = Axis(fig[1, 1]);  ax2 = Axis(fig[2, 1]);
      plot!(ax1, range1, data1; color = :blue, linewidth = 2);
      plot!(ax2, range2, data2; color = :red, linewidth = 2);
      save(figname,fig);
    end

    if sweep_mode == "2D"
      h5name  = "exportFile_"*replace("$(j2)_$(jc1)_$(jc2)","." => "p")*".h5";
      export_to_h5file2D(h5name,data,range1,range2,norm1,norm2);
    elseif sweep_mode == "1D"
      h5name  = "exportFile_"*replace("$(j2)_$(jc1)_$(jc2)","." => "p")*".h5";
      export_to_h5file1D(h5name,data1,data2,range1,range2,norm1,norm2);
    end
    
  end
end

########################################################################################
########################################################################################
########################################################################################