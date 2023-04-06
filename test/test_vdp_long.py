import sys

import pytest

from test.test_vdp import vdp_dormand_prince

@pytest.mark.long
def test_vdp_dormand_prince():
    vdp_dormand_prince(end=101)


# if using 'bazel test ...'
if __name__ == "__main__":
    sys.exit(pytest.main(sys.argv[1:]))
