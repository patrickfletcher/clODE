from math import cos, pi, sqrt
import clode
from typing import List

# Follow clode rhs convention
def sine_curve(
    t: float,
    x_: List[float],
    p_: List[float],
    dx_: List[float],
    aux_: List[float],
    w_: List[float],
) -> None:
    x: float = x_[0]
    dilation: float = p_[0]
    dx: float = cos(t * dilation)
    dx_[0] = dx

# Define the parameters
parameters = {
    "dilation": 1,
}

# Define the initial state
variables = {
    "x": 0,
}

#Activate at t=pi/4, deactivate at t=3pi/2
feature_simulator = clode.FeatureSimulator(
    rhs_equation=sine_curve,
    variables=variables,
    parameters=parameters,
    aux=['dx'],
    observer=clode.Observer.threshold_2,
    stepper=clode.Stepper.rk4,
    dtmax=0.001,
    dt=0.001,
    t_span=(0.0, 4 * pi),
    event_var="x",
    feature_var="x",
    observer_min_x_amp=0.5,
    observer_x_up_thresh=(2+sqrt(2))/4,
    observer_x_down_thresh=0.001,
    observer_dx_down_thresh=0.001,
    observer_dx_up_thresh=0.001,
    observer_max_event_count=100,
    observer_max_event_timestamps=3,
)

# Run the simulation
output = feature_simulator.features()

# Get the number of events (should be 2)
event_count = int(output.get_var_count("event"))
print(event_count)

# Get the timestamps of the events
up_times = output.get_event_data("up")
down_times = output.get_event_data("down")
pass


# Plot the sine curve and visualize the events
import matplotlib.pyplot as plt
import numpy as np

t = np.linspace(0, 4 * pi, 1000)
x = np.sin(t)
plt.plot(t, x)
for event in up_times:
    plt.axvline(x=event, color="green", linestyle="--")

for event in down_times:
    plt.axvline(x=event, color="blue", linestyle="--")
plt.show()
pass
