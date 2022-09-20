import clode

pi = clode.problem_info("samples/lactotroph.cl", 4, 3, 1, 1, ["a", "b", "c", "d"], ["aa", "bb", "cc"], ["dd"])
stepper = "euler"
observer = "basic"
nReps = 1
nPts = 1# 4096
sp = clode.solver_params(0.1, 1.00, 1e-6, 1e-3, 10000000, 10000000, 50)
op = clode.observer_params(0, 0, 100, 1, 1, 0.01, 0.3, 0.2, 0, 0, 1e-7)
open_cl = clode.opencl_resource()
clode_features = clode.clode_features(pi, stepper, observer, True, open_cl, "src/")
tspan = [0.0, 1000.0]
pars = [1.0, 5.0, 1.0] #* nPts
x0 = [0, 0, 0, 0] #* nPts
clode_features.initialize(tspan, x0, pars, sp, op)

clode_features.transient()

print(clode_features.get_tspan())
print(clode_features.get_f())
print(clode_features.get_n_features())
print("Hello world")