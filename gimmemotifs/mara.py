"""Implementation of MARA using pymc."""
import pymc as pm
import numpy as np

def make_model(nmotif, initial_values, real, x):
    """Definition of MARA model."""
    # Prior for motif activities
    a = pm.Normal("a", 0.0, 0.1, size=nmotif, value=initial_values)
    
    @pm.deterministic
    def expected(x=x, a=a):
        """Dot-product of x and a."""
        return np.dot(x, a)
    
    # Prior for sd noise
    sigma = pm.Uniform('sigma', 0.0, 200.0, value=1.0)

    # pylint: disable=unused-variable
    est = pm.Normal('est', mu=expected, tau=sigma**-2, value=real, observed=True)
    
    return locals()
