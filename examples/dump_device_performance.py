''' Benchmark: Time for solution of Lorenz system, 1000 steps of RK4

    For now tests only "transient" with no intermediate storage.
    TODO: test with trajectory and observers
'''

from typing import List
import time
import numpy as np
import matplotlib.pyplot as plt
import clode

def lorenz(
    t: float,
    var: List[float],
    par: List[float],
    dvar: List[float],
    aux: List[float],
    wiener: List[float],
) -> None:
    x: float = var[0]
    y: float = var[1]
    z: float = var[2]

    r: float = par[0]
    s: float = par[1]
    beta: float = par[2]

    dx: float = s*(y - x)
    dy: float = r*x - y - x*z
    dz: float = x*y - beta*z
    dvar[0] = dx
    dvar[1] = dy
    dvar[2] = dz

def test_lorenz_rk4(platform_id, device_id, num_pts, reps):

    parameters = {"r": 27.0, "s": 10.0, "beta": 8.0/3.0}
    variables = {"x": 1.0, "y": 1.0, "z": 1.0}

    t_span = (0.0, 10.0)

    # src_file = "test/lorenz.cl"
    integrator = clode.Simulator(
        # src_file=src_file,
        rhs_equation=lorenz,
        variables=variables,
        parameters=parameters,
        single_precision=True,
        platform_id=platform_id,
        device_id=device_id,
        t_span=t_span,
        stepper=clode.Stepper.rk4,
        dt=0.01,
    )

    t_average = []
    print(f"Lorenz system, 1000 RK4 steps, {reps} repetitions.\nN, min (s), median (s), max (s)")
    for num in num_pts:
        times = []

        integrator.set_repeat_ensemble(num_repeats=num)

        #warm-up passes - seems to make output more reliable
        for _ in range(5):
            integrator.transient(update_x0=False)

        for _ in range(reps):
            t0 = time.perf_counter()
            integrator.transient(update_x0=False)
            times.append(time.perf_counter() - t0)

        # t_mean = sum(times)/reps
        t_min = min(times)
        t_max = max(times)
        t_median = np.median(times)
        print(f"{num}\t {t_min:.3g}, {t_median:.3g}, {t_max:.3g} s")
        t_average.append(t_median)

    return t_average


# if using 'bazel test ...'
if __name__ == "__main__":

    ocl_info = clode.runtime.query_opencl()
    # ocl_info = ocl_info[:3]

    N = [int(2**n) for n in np.arange(0,18)]
    REPS = 50

    device_names = []
    for i, ocl in enumerate(ocl_info):
        if ocl.device_count==0:
            continue
        print("\n", ocl)
        for j, dev in enumerate(ocl.device_info):
            median_time = test_lorenz_rk4(platform_id=i, device_id=j, num_pts=N, reps=REPS)
            plt.plot(N, median_time)
            device_names.append(dev.name)

    plt.legend(device_names)
    plt.xscale('log')
    plt.yscale('log')
    plt.show()
