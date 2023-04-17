import numpy as np

# Define the trapezoid rule for numerical integration
def trapezoid_rule(f, a, b, h):
    n = int((b - a) / h)            # Calculate the number of intervals
    x = np.linspace(a, b, n + 1)    # Create an array of equally spaced points between a and b
    y = f(x)                        # Evaluate the integrand at the x values
    # Calculate the trapezoid rule sum, including the contributions from the first and last points
    return h * (0.5 * y[0] + 0.5 * y[-1] + np.sum(y[1:-1]))

# Define the Riemann sum for numerical integration
def riemann_sum(f, a, b, h):
    n = int((b - a) / h)            # Calculate the number of intervals
    x = np.linspace(a, b, n + 1)    # Create an array of equally spaced points between a and b
    # Evaluate the integrand at the midpoints of each interval (x[:-1] + h / 2)
    y = f(x[:-1] + h / 2)
    # Calculate the Riemann sum by summing the function values and multiplying by the interval width h
    return h * np.sum(y)

# Define Simpson's rule for numerical integration
def simpsons_rule(f, a, b, N):
    if N % 2 != 0:
        raise ValueError("N must be even for Simpson's rule.")
    N = int(N)
    h = (b - a) / N                 # Calculate the interval width
    x = np.linspace(a, b, N + 1)    # Create an array of equally spaced points between a and b
    y = f(x)                        # Evaluate the integrand at the x values
    # Calculate the Simpson's rule sum, including the contributions from the first and last points,
    # the odd-indexed points (multiplied by 4), and the even-indexed points (multiplied by 2)
    return h / 3 * (y[0] + y[-1] + 2 * np.sum(y[2:-1:2]) + 4 * np.sum(y[1:-1:2]))