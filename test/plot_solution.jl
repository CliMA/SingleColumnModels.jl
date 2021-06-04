using Plots
using CSV
# z = CSV.read("test/z1.csv", header=false)[1]
# x = CSV.read("test/x1.csv", header=false)[1]
# x_exact = CSV.read("test/x1_exact.csv", header=false)[1]
# b = CSV.read("test/b1.csv", header=false)[1]

z = CSV.read("test/z2.csv", header=false)[1]
x = CSV.read("test/x2.csv", header=false)[1]
x_exact = CSV.read("test/x2_exact.csv", header=false)[1]
b = CSV.read("test/b2.csv", header=false)[1]
plot(z, x, label="x")
plot!(z, b, label="b")
plot!(z, x_exact, label="x_exact")
savefig("solution.png")
