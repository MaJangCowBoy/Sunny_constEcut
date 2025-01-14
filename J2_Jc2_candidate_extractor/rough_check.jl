include("param.jl");
using CairoMakie

# filename = "LT_minimize_1p0000.dat";
Js_Q = Matrix{Float64}(undef,0,6);

file = open(filename,"r");
for row in eachline(file)
  tmp = parse.(Float64,split(row));
  Js_Q = [Js_Q; tmp[1] tmp[2] tmp[3] tmp[4] tmp[5] tmp[6]];
end;  close(file);

data = zeros(length(J2arr),length(Jc2arr));

for row in eachrow(Js_Q)
  j2 = row[1];  jc1 = row[2];  jc2 = row[3];
  q1 = row[4];  q2  = row[5];  q3  = row[6];
  if     q1 > 0.0 && q2 == 0 # if ordering is q00
    data[findall(J2arr .== j2)[1],findall(Jc2arr .== jc2)[1]] = +q1;
  elseif q1 > 0.0 && q2 > 0 # if ordering is qq0
    data[findall(J2arr .== j2)[1],findall(Jc2arr .== jc2)[1]] = -q1;
  elseif q1 == 0 && q2 == 0 # if ordering is 000
    data[findall(J2arr .== j2)[1],findall(Jc2arr .== jc2)[1]] = NaN;
  end
end

fig = Figure();
axs = Axis(fig[1, 1], xlabel = "J2", ylabel = "Jc2", title = "Q1")
surface!(axs, J2arr, Jc2arr, data[:,:], colormap = :redblue, colorrange = (-1/3,+1/2))
xlims!(axs, 0, 2);  ylims!(axs, -0.5, 0.5);
figname = replace(filename,".dat" => ".png");
save(figname, fig);
