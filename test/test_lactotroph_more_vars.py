# import clode.cpp.clode_cpp_wrapper as clode
# import time
# import numpy as np
#
# pi = clode.problem_info(
#     "test/test.cl", 8, 18, 1, 1, ["V", "n", "m", "b", "h", "h_T", "h_Na", "c"],
#     [
#         'g_CaL', 'g_CaT', 'g_K', 'g_SK', 'g_Kir', 'g_BK', 'g_NaV', 'g_A',
#         'g_leak', 'C_m', 'E_leak', 'tau_m', 'tau_ht', 'tau_n', 'tau_BK',
#         'tau_h', 'tau_hNa', 'k_c'
#     ], ['aux'])
# # pi = clode.problem_info("research/clODE/samples/lactotroph.cl", 4, 3, 1, 1, ["a", "b", "c", "d"], ["aa", "bb", "cc"], ["dd"])
# stepper = "rk4"
# observer = "thresh2"
# nReps = 1
# nPts = 4096
# sp = clode.solver_params(0.5, 1.00, 1e-6, 1e-3, 10000000, 10000000, 50)
# op = clode.observer_params(0, 0, 100, 1, 1, 0.01, 0.3, 0.2, 0, 0, 1e-7)
# open_cl = clode.opencl_resource()
# clode_features = clode.clode_features(pi, stepper, observer, True, open_cl,
#                                       "clode/cpp/")
# tspan = (0.0, 1000.)
# pars = np.array(
#     (1.4, 0, 5, 0, 0, 0, 0, 0, 0.2, 10, -50, 1, 1, 39, 1, 1, 1, 0.03))
#
# pars = np.tile(pars, (nPts, 1)).transpose().flatten()
#
# V = -60.
# n = 0.1
# m = 0.1
# b = 0.1
# h = 0.01
# h_T = 0.01
# h_Na = 0.01
# c = 0.1
#
# x0 = np.array([V, n, m, b, h, h_T, h_Na, c])
# x0 = np.tile(x0, (nPts, 1)).transpose().flatten()
# print(x0)
# print("Init with x0", x0.transpose().flatten().shape)
# print("Init with pars", pars.transpose().flatten().shape)
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
