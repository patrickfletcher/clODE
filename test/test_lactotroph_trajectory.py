# import clode.cpp.clode_cpp_wrapper as clode
# import time
# import numpy as np
#
# # pi = clode.problem_info("test/test.cl", 8, 18, 1, 1,
# #                         ["V", "n", "m", "b", "h", "h_T", "h_Na", "c"],
# #                         ['g_CaL', 'g_CaT', 'g_K', 'g_SK', 'g_Kir', 'g_BK', 'g_NaV', 'g_A', 'g_leak', 'C_m', 'E_leak',
# #                          'tau_m', 'tau_ht', 'tau_n', 'tau_BK', 'tau_h', 'tau_hNa', 'k_c']
# #                         , ['aux'])
# pi = clode.problem_info("samples/lactotroph.cl", 4, 3, 1, 0,
#                         ["v", "n", "f", "c"], ["gcal", "gsk", "gbk"], ["ical"])
# stepper = "dopri5"
# nReps = 1
# nPts = 2
# sp = clode.solver_params(0.5, 100.0, 1e-6, 1e-3, 1000000, 100000, 1)
# open_cl = clode.opencl_resource()
# clode_trajectory = clode.clode_trajectory(pi, stepper, True, open_cl, "clode/cpp/")
# tspan = (0.0, 1000.)
# # pars = np.array((1.4, 0, 5, 0, 0, 0, 0, 0, 0.2, 10, -50, 1, 1, 39, 1, 1, 1, 0.03))
# #
# # pars = np.tile(pars, (nPts, 1)).transpose().flatten()
# #
# # V = -60.
# # n = 0.1
# # m = 0.1
# # b = 0.1
# # h = 0.01
# # h_T = 0.01
# # h_Na = 0.01
# # c = 0.1
# #
# # x0 = np.array([V, n, m, b, h, h_T, h_Na, c])
# # x0 = np.tile(x0, (nPts, 1)).transpose().flatten()
# pars = [1.5] * nPts + [3.0] * nPts + [1.0] * nPts
# x0 = [0, 0, 0, 0] * nPts
# clode_trajectory.build_cl()
# clode_trajectory.initialize(tspan, x0, pars, sp)
#
# mySeed = 1
# clode_trajectory.seed_rng(mySeed)
# clode_trajectory.transient()
# # clode_trajectory.shift_x0()
#
# start_time = time.perf_counter_ns()
#
# for _ in range(nReps):
#     clode_trajectory.trajectory()
#
# end_time = time.perf_counter_ns()
#
# n_stored = clode_trajectory.get_n_stored()
# tt = clode_trajectory.get_t()
# xx = clode_trajectory.get_x()
# xf = clode_trajectory.get_x0()
# trajIx = 0
#
# print(f"\ntf , F:")
# for ii in range(n_stored[trajIx]):
#     print(f"{tt[ii * nPts + trajIx]:.3f}", end=" ")
#     for jj in range(4):
#         print(f"{xx[ii * nPts * 4 + jj * nPts + trajIx]:.3f}", end=" ")
#
#     print()
#
# # std::cout<< "Timepoints stored: " << nStored[trajIx] << "/" << nStoreMax << "\n";
# print(f"Timepoints stored: {n_stored[trajIx]}/{max(n_stored)}")
#
# elapsed_time = round((end_time - start_time) / 10**6, 3)
# print(f"Compute time: {elapsed_time}ms")
