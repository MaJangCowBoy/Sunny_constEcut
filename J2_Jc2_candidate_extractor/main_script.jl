#? Please execute with "julia_wh --threads=auto"
#? or "julia_wh --threads=4" etc.

using Statistics, LinearAlgebra, Printf
include("param.jl");  include("function_bundle.jl");

npar = Threads.nthreads();

##################################################################

Js_arr = Matrix{Float64}(undef,0,4);
Qs_arr = Matrix{Float64}(undef,0,3);

# Jc1Jc2_rows = [1.20 0.40; 1.20 0.50; 1.20 0.60];
# for (Jc1,Jc2) in eachrow(Jc1Jc2_rows)
## If you want a multiple results for (Jc1,Jc2) combinations, 
## then turn on

for j2 in J2arr, jc2 in Jc2arr, jc1 in Jc1
  F = jc1 + jc2;  G = jc1 - jc2 * 0.5;
  j3 = 1/2 * (1 - (F*G - jc2*F/2 + 2*jc2*G )/√(F*F-2*F*G+4*G*G));
  if isnan(j3)  j3 = 0.5  end
  global Js_arr = [Js_arr; j2 j3 jc1 jc2];
end

q00_grid = range(0,1/2,501);
qq0_grid = range(0,1/3,334);  Q3 = 0.0;
for q in q00_grid  global Qs_arr = [Qs_arr; q 0 Q3];  end
for q in qq0_grid  global Qs_arr = [Qs_arr; q q Q3];  end

LenJs = size(Js_arr,1);  Js_Q  = zeros(LenJs,6);
LenQs = size(Qs_arr,1);

Threads.@threads for id in 1:npar

  SttPt = 1 + (id-1) * div(LenJs,npar);
  EndPt = id * div(LenJs,npar);
  if id == npar  EndPt = LenJs;  end

  for it1 = SttPt:EndPt
    tmp = zeros(LenQs);
    for it2 = 1:LenQs
      tmp[it2] = calc_matrix(Js_arr[it1,:],Qs_arr[it2,:]);
    end;  MinPos = find_minpos(tmp,10*eps());
    
    Qs_tobe = zeros(length(MinPos),3);
    for it2 in eachindex(MinPos)
      Qs_tobe[it2,:] = Qs_arr[MinPos[it2],:];
    end

    Qs_sort = sortslices(Qs_tobe,dims=1);
    Qs = Qs_sort[1,:];  
    Q1 = Qs[1];  Q2 = Qs[2];  Q3 = Qs[3];
    J2  = Js_arr[it1,1];  J3  = Js_arr[it1,2];
    Jc1 = Js_arr[it1,3];  Jc2 = Js_arr[it1,4];

    Js_Q[it1,:] = [J2,Jc1,Jc2,Q1,Q2,Q3];
    # Js_Q[it1,:] = [Js_arr[it1,:],Qs];
  end
end

file = open(filename,"w");
for rows in eachrow(Js_Q)
  for elem in rows  print(file, @sprintf("%.6f ",elem))  end
  println(file);
end;  close(file);

# end
## If you want a multiple results for (Jc1,Jc2) combinations,
## then turn on

##################################################################