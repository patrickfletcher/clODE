import pytest

import clode


@pytest.mark.parametrize(
    "device_type, vendor, platform_id, device_id, device_ids",
    [
        [
            clode.CLDeviceType.DEVICE_TYPE_CPU,
            clode.CLVendor.VENDOR_ANY,
            0,
            None,
            None,
        ],
        [
            clode.CLDeviceType.DEVICE_TYPE_CPU,
            clode.CLVendor.VENDOR_ANY,
            None,
            0,
            None,
        ],
        [
            clode.CLDeviceType.DEVICE_TYPE_CPU,
            clode.CLVendor.VENDOR_ANY,
            None,
            None,
            [0],
        ],
        [
            clode.CLDeviceType.DEVICE_TYPE_CPU,
            None,
            0,
            None,
            None,
        ],
        [
            clode.CLDeviceType.DEVICE_TYPE_CPU,
            None,
            None,
            0,
            None,
        ],
        [
            clode.CLDeviceType.DEVICE_TYPE_CPU,
            None,
            None,
            None,
            [0],
        ],
        [
            None,
            clode.CLVendor.VENDOR_ANY,
            0,
            None,
            None,
        ],
        [
            None,
            clode.CLVendor.VENDOR_ANY,
            None,
            0,
            None,
        ],
        [
            None,
            clode.CLVendor.VENDOR_ANY,
            None,
            None,
            [0],
        ],
        [
            None,
            None,
            0,
            0,
            [0],
        ],
    ],
)
def test_init_features_runtime_with_incorrect_config_fails(
    device_type, vendor, platform_id, device_id, device_ids
):
    input_file: str = "test/van_der_pol_oscillator.cl"

    t_span = (0.0, 1000.0)

    with pytest.raises(ValueError):
        _ = clode.TrajectorySimulator(
            src_file=input_file,
            variables={"x": 1.0, "y": 1.0},
            parameters={"mu": 1.0},
            num_noise=0,
            stepper=clode.Stepper.dormand_prince,
            t_span=t_span,
            device_type=device_type,
            vendor=vendor,
            platform_id=platform_id,
            device_id=device_id,
            device_ids=device_ids,
        )

    with pytest.raises(ValueError):
        _ = clode.FeatureSimulator(
            src_file=input_file,
            variables={"x": 1.0, "y": 1.0},
            parameters={"mu": 1.0},
            num_noise=0,
            stepper=clode.Stepper.dormand_prince,
            t_span=t_span,
            device_type=device_type,
            vendor=vendor,
            platform_id=platform_id,
            device_id=device_id,
            device_ids=device_ids,
        )


def test_fast_and_slow():
    from typing import List

    import matplotlib.pyplot as plt
    import numpy as np

    import clode

    # The FitzHugh-Nagumo model

    def fitzhugh_nagumo(
        time: float,
        variables: List[float],
        a: float,
        b: float,
        current: float,
        epsilon: float,
    ) -> List[float]:
        v: float = variables[0]
        w: float = variables[1]
        dv = v - v**3 / 3 - w + current
        dw = epsilon * (v + a - b * w)
        return [dv, dw]

    a = [0.7, 0.8, 0.9, 1.0]

    simulator = clode.FeatureSimulator(
        rhs_equation=fitzhugh_nagumo,
        variables={"v": np.arange(-2, 2, 0.2), "w": -1.0},
        parameters={"a": a, "b": 0.8, "epsilon": 0.01, "current": 1.0},
    )
    print("Foo")
