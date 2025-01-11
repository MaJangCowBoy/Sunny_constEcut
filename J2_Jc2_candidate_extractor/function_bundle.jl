include("param.jl");

##################################################################

function calc_matrix(J2J3JcJd, Q1Q2Q3)
  J1 = 1.0;
  J2  = J2J3JcJd[1];  J3  = J2J3JcJd[2];  
  Jc1 = J2J3JcJd[3];  Jc2 = J2J3JcJd[4];
  Q1 = Q1Q2Q3[1];  Q2 = Q1Q2Q3[2];  Q3 = Q1Q2Q3[3];
  Q = hcat(Q1,Q2,Q3);

  Emat = zeros(ComplexF64,2,2);  
  Emat[1,1] =  J1 * sum( exp.(2π*im*Q* dAA_1NN) ) + 
               J2 * sum( exp.(2π*im*Q* dAA_2NN) ) +
               J3 * sum( exp.(2π*im*Q* dAA_3NN) ) ;
  Emat[1,2] = Jc1 * sum( exp.(2π*im*Q*dAB_c1NN) ) + 
              Jc2 * sum( exp.(2π*im*Q*dAB_c2NN) ) ;
  Emat[2,1] = Jc1 * sum( exp.(2π*im*Q*dBA_c1NN) ) +
              Jc2 * sum( exp.(2π*im*Q*dBA_c2NN) ) ;
  Emat[2,2] =  J1 * sum( exp.(2π*im*Q* dBB_1NN) ) + 
               J2 * sum( exp.(2π*im*Q* dBB_2NN) ) +
               J3 * sum( exp.(2π*im*Q* dBB_3NN) ) ;
  k = eigen(Emat);
  E1 = real(k.values[1]);  E2 = real(k.values[2]);

  return minimum([E1,E2]);
end

##################################################################

function find_minpos(arr, tolerance)
  min_val = minimum(arr);  
  min_pos = [];
  for (i, val) in enumerate(arr)
    if abs(val - min_val) <= tolerance  push!(min_pos, i)  end
  end;  min_pos = convert(Array{Int64},min_pos);

  return min_pos;
end

##################################################################