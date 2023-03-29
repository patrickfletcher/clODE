# import clode.cpp.clode_cpp_wrapper as clode
# import time
#
# pi = clode.problem_info("samples/lactotroph.cl", 4, 3, 1, 1,
#                         ["a", "b", "c", "d"], ["aa", "bb", "cc"], ["dd"])
# #pi = clode.problem_info("research/clODE/samples/lactotroph.cl", 4, 3, 1, 1, ["a", "b", "c", "d"], ["aa", "bb", "cc"], ["dd"])
# stepper = "euler"
# observer = "thresh2"
# nReps = 1
# nPts = 4096
# sp = clode.solver_params(0.1, 1.00, 1e-6, 1e-3, 10000000, 10000000, 50)
# op = clode.observer_params(0, 0, 100, 1, 1, 0.01, 0.3, 0.2, 0, 0, 1e-7)
# open_cl = clode.opencl_resource()
# clode_features = clode.clode_features(pi, stepper, observer, True, open_cl,
#                                       "clode/cpp/")
# tspan = [0.0, 1000.0]
# pars = [1.0, 5.0, 1.0] * nPts
# x0 = [0, 0, 0, 0] * nPts
# clode_features.build_cl()
# clode_features.initialize(tspan, x0, pars, sp, op)
#
# mySeed = 1
# clode_features.seed_rng(mySeed)
# clode_features.transient()
# clode_features.shift_x0()
#
# start_time = time.perf_counter_ns()
#
# for _ in range(nReps):
#     clode_features.features()
#
# end_time = time.perf_counter_ns()
#
# F = clode_features.get_f()
# tspan = clode_features.get_tspan()
# print(clode_features.get_feature_names(), clode_features.get_n_features())
# print(f"\ntf={tspan[0]} , F:")
# for i in range(clode_features.get_n_features()):
#     print(F[i * nPts], end=" ")
#
# print()
#
# elapsed_time = round((end_time - start_time) / 10**6, 3)
# print(f"Compute time: {elapsed_time}ms")
