using DelimitedFiles, CairoMakie, Printf
include("param.jl");

# filename = "LT_minimize_0p0000.dat";
data = Matrix{Float64}(undef, 0, 6);

f = open(filename, "r")
for line in eachline(f)
  tok = parse.(Float64,split(line))
  global data = [data; tok']
end;  close(f);

phasediagram = NaN * zeros(length(J2arr),length(Jc2arr));

for row in eachrow(data)
  if row[4] > 0 && abs(row[5]) < 1e-5
    x = findall(J2arr .== row[1])[1];
    global jc1 = row[2];
    y = findall(Jc2arr .== row[3])[1];
    phasediagram[x,y] = 1;
  end
end

fig = Figure();
ax = Axis(fig[1, 1]);
heatmap!(ax, J2arr, Jc2arr, phasediagram);
xlims!(ax, minimum(J2arr), maximum(J2arr));  ylims!(ax, minimum(Jc2arr), maximum(Jc2arr));
figname = replace(filename, ".dat" => ".png");
save(figname, fig);

phasediagram_trim = copy(phasediagram);
for id1 in axes(phasediagram,1), id2 in axes(phasediagram,2)
  if isnan(phasediagram[id1,id2])
    for x = (id1-1):(id1+1), y = (id2-1):(id2+1)
      if x > 0 && x < size(phasediagram,1)+1 && y > 0 && y < size(phasediagram,2)+1
        phasediagram_trim[x,y] = NaN;
      end
    end
  end
end

fig = Figure();
ax = Axis(fig[1, 1]);
heatmap!(ax, J2arr, Jc2arr, phasediagram_trim);
xlims!(ax, minimum(J2arr), maximum(J2arr));  ylims!(ax, minimum(Jc2arr), maximum(Jc2arr));
figname = replace(filename, ".dat" => "_trim.png");
save(figname, fig);

filename = replace(filename, ".dat" => "_trim.dat");
f = open(filename, "w")
for id1 in axes(phasediagram_trim,1), id2 in axes(phasediagram_trim,2)
  if !isnan(phasediagram_trim[id1,id2])
    j2 = J2arr[id1];  jc2 = Jc2arr[id2];
    str = @sprintf("%.3f %.3f %.3f %.3f %.3f %.3f", j2, jc1, jc2, 0.333, 0.0, 0.0);
    println(f, str)
  end
end;  close(f);

