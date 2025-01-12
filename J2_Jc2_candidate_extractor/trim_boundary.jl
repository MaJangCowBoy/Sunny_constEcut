using DelimitedFiles, CairoMakie, Printf

data = Matrix{Float64}(undef, 0, 6);
f = open("LT_minimize_0p0000.dat", "r")
for line in eachline(f)
  tok = parse.(Float64,split(line))
  data = [data; tok']
end;  close(f);

phasediagram = NaN * zeros(101,51);

for row in eachrow(data)
  if row[4] > 0 && abs(row[5]) < 1e-5
    x = Int(round( 1 + row[1]/0.02));
    y = Int(round(26 + row[3]/0.02));
    phasediagram[x,y] = 1;
  end
end

fig = Figure(resolution = (800, 400));
ax = Axis(fig[1, 1]);
heatmap!(ax, 0.0:0.02:2.0, 0-0.5:0.02:0.5, phasediagram);
xlims!(ax, 0.0, 2.0);  ylims!(ax, -0.5, 0.5);
save("LT_minimize_0p0000.png", fig);

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

fig = Figure(resolution = (800, 400));
ax = Axis(fig[1, 1]);
heatmap!(ax, 0.0:0.02:2.0, 0-0.5:0.02:0.5, phasediagram_trim);
xlims!(ax, 0.0, 2.0);  ylims!(ax, -0.5, 0.5);
save("LT_minimize_0p0000_trim.png", fig);

f = open("LT_minimize_0p0000_trim.dat", "w")
for id1 in axes(phasediagram_trim,1), id2 in axes(phasediagram_trim,2)
  if !isnan(phasediagram_trim[id1,id2])
    j2 = (id1-1)*0.02;  jc2 = (id2-26)*0.02;
    str = @sprintf("%.3f %.3f %.3f %.3f %.3f %.3f", j2, 0.000, jc2, 0.333, 0.0, 0.0);
    println(f, str)
  end
end;  close(f);

