using DataFrames, CSV, StatsPlots


f1=DataFrame(CSV.File("lapacke_single_precision.dat", delim=" ", header=false))
f2=DataFrame(CSV.File("lapacke_double_precision.dat", delim=" ", header=false))
f3=DataFrame(CSV.File("tuvx_single_precision.dat", delim=" ", header=false))
f4=DataFrame(CSV.File("tuvx_double_precision.dat", delim=" ", header=false))


lapack_errors_single_precision = f1[:, 1]
lapack_times_single_precision = f1[:, 2]

tuvx_errors_single_precision = f3[:, 1]
tuvx_times_single_precision = f3[:, 2]

lapack_errors_double_precision = f2[:, 1]
lapack_times_double_precision = f2[:, 2]

tuvx_errors_double_precision = f4[:, 1]
tuvx_times_double_precision = f4[:, 2]

nam = repeat([500, 1000, 10000, 100000, 1000000], outer=2)
# julia has a weird way of ordering groups in bar plots (need to make it better)
nam = [ n > 500 ? "$n" : 
        "  $n" 
        for n in nam]

data = [lapack_errors_single_precision; tuvx_errors_single_precision]
groups =  repeat(["LAPACKE", "TUV-X"], inner = 5)
plot(groupedbar(nam, data, groups=groups))
title!("Comparing Accuracy (Single Precision)")
xlabel!("System Size")
ylabel!("Relative Error")
savefig("single_precision_errors.png")

data = [lapack_times_single_precision; tuvx_times_single_precision]
groups =  repeat(["LAPACKE", "TUV-X"], inner = 5)
plot(groupedbar(nam, data, groups=groups))
title!("Comparing Speed (Double Precision)")
xlabel!("System Size")
ylabel!("Run Time (ms)")
savefig("single_precision_times.png")

data = [lapack_errors_double_precision; tuvx_errors_double_precision]
groups =  repeat(["LAPACKE", "TUV-X"], inner = 5)
plot(groupedbar(nam, data, groups=groups))
title!("Comparing Accuracy (Double Precision)")
xlabel!("System Size")
ylabel!("Relative Error")
savefig("double_precision_errors.png")

data = [lapack_times_double_precision; tuvx_times_double_precision]
groups =  repeat(["LAPACKE", "TUV-X"], inner = 5)
plot(groupedbar(nam, data, groups=groups))
title!("Comparing Speed (Double Precision)")
xlabel!("System Size")
ylabel!("Run Time (ms)")
savefig("double_precision_times.png")

