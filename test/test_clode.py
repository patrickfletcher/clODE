import clode
pi = clode.problem_info("samples/lactotroph.cl", 4, 3, 1, 1, ["a", "b", "c", "d"], ["aa", "bb", "cc"], ["dd"])


stepper = "euler"
observer = "basic"

#//parameters for solver and objective function
#std::vector<double> tspan({0.0,1000.0});

nReps = 1
nPts = 1# 4096

sp = clode.solver_params(0.1, 1.00, 1e-6, 1e-3, 10000000, 10000000, 50)

# SolverParams<double> sp;
# sp.dt=0.1;
# sp.dtmax=1.00;
# sp.abstol=1e-6;
# sp.reltol=1e-3;
# sp.max_steps=10000000;
# sp.max_store=10000000;
# sp.nout=50;

op = clode.observer_params(0, 0, 100, 1, 1, 0.01, 0.3, 0.2, 0, 0, 1e-7)
# op.eVarIx=0;
# op.fVarIx=0;
# op.maxEventCount=100;
# op.minXamp=1;
# op.nHoodRadius=0.01;
# op.xUpThresh=0.3;
# op.xDownThresh=0.2;
# op.dxUpThresh=0;
# op.dxDownThresh=0;
# op.eps_dx=1e-7;

open_cl = clode.opencl_resource()

clode_features = clode.clode_features(pi, stepper, observer, True, open_cl, "src/")

tspan = [0.0, 1000.0]

pars = [1.0, 5.0, 1.0] #* nPts
x0 = [0, 0, 0] #* nPts

clode_features.initialize(tspan, x0, pars, sp, op)

print("Hello world")