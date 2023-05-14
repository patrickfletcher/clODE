import sys
from test.test_vdp import vdp_dormand_prince

import pytest

from clode import convert_xpp_file


def test_vdp_converter_converts_vdp():
    input_file = "test/xpp/van_der_pol_oscillator.xpp"
    reference_file = "test/xpp/van_der_pol_oscillator_reference.cl"

    output_file = convert_xpp_file(input_file)

    output_content = open(output_file, "r").read()
    reference_content = open(reference_file, "r").read()

    assert output_content == reference_content

    vdp_dormand_prince(end=5, input_file=input_file)


# if using 'bazel test ...'
if __name__ == "__main__":
    sys.exit(pytest.main(sys.argv[1:]))
