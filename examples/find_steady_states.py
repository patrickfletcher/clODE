import clode
import numpy as np

# exponential decay to a specified value
def decay(t:float, y:list[float], p:list[float], dy:list[float], aux:list[float], w:list[float]) -> None:
    y_target:float = p[0]
    tau:float = p[1]
    dy[0] = (y_target - y[0])/tau

variables = {"y": 0.0}
parameters = {"y_target": 5.0, "tau": 20.0}

simulator = clode.Simulator(rhs_equation=decay, variables=variables, parameters=parameters)

num_sims = 6
tau_values = 10. + 90.*np.random.random(size=(num_sims,))
y0 = 10.*np.random.random(size=(num_sims,))
simulator.set_ensemble(parameters={"tau":tau_values}, variables={"y":y0})

print(f"tau: {tau_values.transpose()}")
print(f"y0:  {y0.transpose()}")

t_span = (0., 100.)
maxdiff = 5.
while maxdiff>1e-3:
    ylast = simulator.get_initial_state()
    yf = simulator.transient(t_span = t_span, update_x0=True, fetch_results=True)
    maxdiff = np.max(np.abs(yf - ylast))
    
    print(f"yf:  {yf[:,0].transpose()}")
    # print(f"{maxdiff=}")