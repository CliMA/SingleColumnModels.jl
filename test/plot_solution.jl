using Plots
using CSV
z = CSV.read("test/z.csv", header=false)[1]
x = CSV.read("test/x.csv", header=false)[1]
x_exact = CSV.read("test/x_exact.csv", header=false)[1]
b = CSV.read("test/b.csv", header=false)[1]
plot(z, x, label="x")
plot!(z, b, label="b")
plot!(z, x_exact, label="x_exact")
savefig("solution.png")
