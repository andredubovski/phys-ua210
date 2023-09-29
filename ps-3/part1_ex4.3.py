import numpy as np

def my_function(input_value=None, data_type=float):
    return input_value * (input_value - 1)

def my_analytical_derivative(input_value=None, data_type=float):
    return 2 * input_value - 1

def my_numerical_derivative(input_value=None, step_size=None, data_type=float):
    return (my_function((input_value + step_size)) - my_function(input_value)) / step_size

print(my_analytical_derivative(1))

for i in range(1,8):
    step_size = 10**(-2 * i)
    print(my_numerical_derivative(1, step_size))
