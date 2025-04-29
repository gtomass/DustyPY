import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import numpy as np
from DustyPY.MCfit import fit
from DustyPY.Data import Data
import matplotlib.pyplot as plt

# Define the exponential model as a top-level function
def exp_model(theta, data):
    """
    Exponential model: y = theta[0] * exp(-theta[1] * x)
    """
    x = data.xdata[0]
    return theta[0] * np.exp(-theta[1] * x)

# Define the sum-of-squares function as a top-level function
def sos_function(theta, data):
    """
    Sum-of-squares function for the exponential model.
    """
    ydata = data.ydata[0]
    ymodel = exp_model(theta, data)
    return np.sum((ydata - ymodel) ** 2)

def test_fit_initialization():
    """
    Test the initialization of the fit class.
    """
    # Create synthetic data
    xdata = np.linspace(0, 10, 50)
    ydata = 3.0 * np.exp(-0.5 * xdata)
    data = Data(xdata=xdata, ydata=ydata)

    # Initialize the fit object
    mcmc_fit = fit(data=data)

    assert mcmc_fit.get_Data() == data, "Data not set correctly"

def test_fit_set_parameters():
    """
    Test setting parameters in the fit class.
    """
    # Create synthetic data
    xdata = np.linspace(0, 10, 50)
    ydata = 3.0 * np.exp(-0.5 * xdata)
    data = Data(xdata=xdata, ydata=ydata)

    # Initialize the fit object
    mcmc_fit = fit(data=data)

    # Set parameters
    param = {
        'theta1': {'theta0': 3.0, 'minimum': 0.0, 'maximum': 10.0, 'sample': True},
        'theta2': {'theta0': 0.5, 'minimum': 0.0, 'maximum': 2.0, 'sample': True},
    }
    mcmc_fit.set_Param(param)

    assert mcmc_fit.get_Param() == param, "Parameters not set correctly"

def test_fit_run_simulation():
    """
    Test running the MCMC simulation in the fit class.
    """
    # Create synthetic data
    xdata = np.linspace(0, 10, 50)
    ydata = 3.0 * np.exp(-0.5 * xdata) + np.random.normal(0, 0.1, len(xdata))
    data = Data(xdata=xdata, ydata=ydata)

    # Initialize the fit object
    mcmc_fit = fit(data=data, ncpus=4)

    # Set parameters
    param = {
        'theta1': {'theta0': 3.0, 'minimum': 0.0, 'maximum': 10.0, 'sample': True},
        'theta2': {'theta0': 0.5, 'minimum': 0.0, 'maximum': 2.0, 'sample': True},
    }
    mcmc_fit.set_Param(param)

    # Run the MCMC simulation
    mcmc_fit.fit(Chi2=sos_function)

    # Retrieve results
    results = mcmc_fit.get_Results()
    stats = mcmc_fit.get_Stats()

    # Assertions to verify results
    assert results is not None, "Results should not be None"
    assert stats is not None, "Stats should not be None"
    assert 'chain' in results, "Results should contain the chain"
    assert results['chain'].shape[0] > 0, "Chain should have samples"
    assert results['chain'].shape[1] == len(param), "Chain should have the correct number of parameters"

def test_fit_plot_results():
    """
    Test plotting results in the fit class.
    """
    # Create synthetic data
    xdata = np.linspace(0, 10, 50)
    ydata = 3.0 * np.exp(-0.5 * xdata) + np.random.normal(0, 0.1, len(xdata))
    data = Data(xdata=xdata, ydata=ydata)

    # Initialize the fit object
    mcmc_fit = fit(data=data)

    # Set parameters
    param = {
        'theta1': {'theta0': 3.0, 'minimum': 0.0, 'maximum': 10.0, 'sample': True},
        'theta2': {'theta0': 0.5, 'minimum': 0.0, 'maximum': 2.0, 'sample': True},
    }
    mcmc_fit.set_Param(param)

    # Run the MCMC simulation
    mcmc_fit.fit(Chi2=sos_function)

    # Plot results
    chain = mcmc_fit.get_Results()['chain']
    burnin = int(mcmc_fit.get_Results()['nsimu'] / 2)
    plt.figure(figsize=(10, 5))
    for i, name in enumerate(mcmc_fit.get_Results()['names']):
        plt.plot(chain[burnin:, i], label=name)
    plt.xlabel("Iteration")
    plt.ylabel("Parameter Value")
    plt.legend()
    plt.title("MCMC Chains")
    plt.show()

if __name__ == "__main__":
    test_fit_initialization()
    test_fit_set_parameters()
    test_fit_run_simulation()
    test_fit_plot_results()
