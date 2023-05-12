from E2.E2 import trapezoid_rule, riemann_sum, simpsons_rule
import numpy as np
##################################################
# Importing E2.py as a module
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

##################################################

# Define the integrands for I1 and I2


def integrand_I1(x):
    return np.exp(-x) / x


def integrand_I2(x):
    # return x * np.sin(1 / x)   This line would return a division by zero error
    # To mitigate this, we can use a boolean mask
    result = np.zeros_like(x)
    mask = x != 0
    result[mask] = x[mask] * np.sin(1 / x[mask])
    return result


# Initialize the parameters for the integrals
a1, b1 = 1, 100
a2, b2 = 0, 1
h_initial_1 = b1 - a1
h_initial_2 = b2 - a2

# Function to calculate the integral using different methods and halving the interval width
# Note that Simpson's rule has to be treated differently


def calculate_integral(integrand, a, b, initial_h, integration_method, target_relative_change=1e-4):
    h = initial_h if integration_method is not simpsons_rule else 2
    previous_result = integration_method(integrand, a, b, h)
    while True:
        h = h/2 if integration_method is not simpsons_rule else h*2
        current_result = integration_method(integrand, a, b, h)
        relative_change = np.abs(
            (current_result - previous_result) / previous_result)
        if relative_change < target_relative_change:
            break
        previous_result = current_result
    return (current_result, h)


# Calculate the integrals I1 and I2 with different methods
for method_name, method in [("Trapezoid rule", trapezoid_rule), ("Riemann sum", riemann_sum), ("Simpson's rule", simpsons_rule)]:
    result_I1 = calculate_integral(integrand_I1, a1, b1, h_initial_1, method)
    result_I2 = calculate_integral(integrand_I2, a2, b2, h_initial_2, method)
    print(f"{method_name}:")
    print(
        f"I1 ≈ {result_I1[0]:.6f}, at {'N' if method is simpsons_rule else 'h'} = {format(result_I1[1], '.4g')}")
    print(
        f"I2 ≈ {result_I2[0]:.6f}, at {'N' if method is simpsons_rule else 'h'} = {format(result_I2[1], '.4g')}\n")
    # format(number, '.4g') rounds a number to 4 significant digits
