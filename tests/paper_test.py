"""
This test is made to be run using pytest to verify, that the results of the v1 paper are reproduceable
using the current version of quadtune.

This test will not work if it is not executed either from the quadtest directory or using pytest
"""

import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from tuning_files import quadtune_driver
import numpy as np

def test_paper_tuning():
    
    parameter_values = quadtune_driver.main(["-c", "tuning_files.configs_tests.config_v1_paper_test"])
    true_parameter_values = [0.7970267, 0.4655446, 0.0880349, 0.0004552707, 0.1357132]
    assert np.allclose(parameter_values.flatten(), true_parameter_values, rtol=1e-5)

if __name__ == "__main__":
    test_paper_tuning()