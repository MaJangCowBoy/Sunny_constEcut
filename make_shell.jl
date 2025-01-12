using Printf

data = Matrix{Float64}(undef, 0, 6);
f = open("LT_minimize_0p0000_trim.dat", "r");
for line in eachline(f)
  tok = parse.(Float64,split(line))
  data = [data; tok']
end;  close(f);


f = open("autoRun.sh","w");

println(f, "#!/bin/bash")
for (id,row) in enumerate(eachrow(data))
  j2 = row[1];  jc2 = row[3];
  str = @sprintf("julia_wh main_script.jl 3Q 1D %.3f %.3f", j2, jc2);
  println(f, str)
  if id%104 == 0
    println(f, "wait")
  end
end;
println(f, "wait")
println(f, "echo 'All jobs are done.'")
println(f, "exit 0")

close(f);