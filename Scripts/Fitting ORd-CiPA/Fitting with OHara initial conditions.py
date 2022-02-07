# -*- coding: utf-8 -*-
"""~

@author: yanral
"""
import sabs_pkpd
import pints
import numpy as np


# Load the MMT model for the O'Hara CiPA model
filename = '/pstore/home/barraly/Codes/Gamma_0 impact on fitting/Ohara CiPA - algebraic voltage.mmt'
s = sabs_pkpd.load_model.load_simulation_from_mmt(filename)
s.set_tolerance(1e-08, 1e-08)
default_state = s.state()

# Save the initial conditions published in OHara CiPA model
Ohara_init_conds = default_state.copy()


# Define the functions to make sure there is consistency between the initial conditions
def G0_calc(Ki = 144.65559, Kss = 144.65556, Nai = 7.268, Nass = 7.26809,
            Cai = 8.6e-5, Cansr = 1.61957, Cajsr = 1.571234014, Cass = 8.49e-5,
            V=-88, extraK = 5.4, extraNa = 140, extraCa = 1.8):
    tot_cai = Cai * (1 + 0.05 / (Cai + 0.00238) + 0.07/(Cai + 0.0005))
    tot_cass = Cass * (1 + 0.047 / (Cass + 0.00087) + 1.124/(Cass + 0.0087))
    tot_cajsr = Cajsr * (1 + 10 / (Cajsr + 0.8))
    return -V / (96485 * 2.583592e-05) * 0.0001533576 + Ki + Kss * 0.029411764705882353 + Nai + Nass * 0.029411764705882353 + 2*(tot_cai + tot_cass * 0.029411764705882353 + Cansr * 0.08117647059 + tot_cajsr * 0.007059) - extraK - extraNa - 2 * extraCa

def Ki_calc(G0, Nai = 7.268, Nass = 7.26809, Cai = 8.6e-5, Cansr = 1.61957, Cajsr = 1.571234014, Cass = 8.49e-5, V=-88, extraK = 5.4, extraNa = 140, extraCa = 1.8):
    tot_cai = Cai * (1 + 0.05 / (Cai + 0.00238) + 0.07/(Cai + 0.0005))
    tot_cass = Cass * (1 + 0.047 / (Cass + 0.00087) + 1.124/(Cass + 0.0087))
    tot_cajsr = Cajsr * (1 + 10 / (Cajsr + 0.8))
    return (V / (96485 * 2.583592e-05) * 0.0001533576 + extraK + extraNa + 2 * extraCa + G0 - Nai - Nass * 0.029411764705882353 - 2*(tot_cai + tot_cass * 0.029411764705882353 + Cansr * 0.08117647059 + tot_cajsr * 0.007059)) / 1.029411764705882353 

def compute(Gamma_0):
    # Reinitialise the myokit.Simulation
    s.reset()
    
    # Set the initial conditions for Ki and Kss so that the initial conditions match with the value of Gamma_0
    initial_state = default_state.copy()
    initial_K = Ki_calc(Gamma_0,
                       Nai = default_state[1],
                       Nass = default_state[2],
                       Cai = default_state[5],
                       Cansr = default_state[7],
                       Cajsr = default_state[8],
                       Cass = default_state[6])
    initial_state[3] = initial_K
    initial_state[4] = initial_K
    s.set_state(initial_state)
    
    # Record the action potential at the limit cycle
    s.pre(2000000)
    out = s.run(1000, log_interval = 1)
    
    return out['membrane.V']


# Run the model with the published original initial conditions and the Gamma_0 value associated with it
Original_Gamma_0 = 7.801116
data_to_fit = compute(Original_Gamma_0)


# Define the time points on which to read the voltage
time_points = np.linspace(0, 999, 1000)

# Define the fitted parameters and initial point
parameters_to_fit = ['ical.rescale', 'ikr.rescale', 'IKs.rescale', 'INa.rescale', 'INaL.rescale']
true_values = np.array([1, 1, 1, 1, 1])
x0 = np.random.uniform(low = 0.2, high = 5, size = 5)
print('Initial point for fitting : ' + str(x0))


# Set the boundaries for fitting
mini = np.ones(5) * 0.1
print('Lower boundary : ' + str(mini))
maxi = np.ones(5) * 10
print('Upper boundary : ' + str(maxi))


# Prepare Pints optimisation routine
# In[Run the optimisation]
# Run optimisation of the 13 params
class MyModel(pints.ForwardModel):
    def n_parameters(self):
        # Define the amount of fitted parameters
        return sabs_pkpd.constants.n

    def simulate(self, parameters, times):
        sabs_pkpd.constants.n = len(parameters)
        
        # Set the rescaling parameters
        for p, label in enumerate(parameters_to_fit):
            s.set_constant(label, parameters[p])
            
        # In case there is a numerical error
        try:
            out = compute(Gamma_0 = Original_Gamma_0)
            return out
        
        except:
            print('Simulation error. Continuing anyway...')
            return np.zeros(len(times))

# Define the parameters and method for fitting
sigma0 = [0.1, 0.1, 0.1, 0.1, 0.1]
sabs_pkpd.constants.n = len(parameters_to_fit)


problem = pints.SingleOutputProblem(
    model=MyModel(),
    times=time_points,
    values=data_to_fit)
boundaries = pints.RectangularBoundaries(mini, maxi)
function = pints.SumOfSquaresError(problem)
method = pints.CMAES


# Detailed code from PINTS controller to log the CSV files during fitting
# Convert x0 to vector
# This converts e.g. (1, 7) shapes to (7, ), giving users a bit more
# freedom with the exact shape passed in. For example, to allow the
# output of LogPrior.sample(1) to be passed in.
x0 = pints.vector(x0)

# Check if minimising or maximising
minimising = not isinstance(function, pints.LogPDF)

# Store transform for later detransformation: if using a transform, any
# parameters logged to the filesystem or printed to screen should be
# detransformed first!
transform = None

# Create optimiser
optimiser = method(x0, sigma0, boundaries)

# Check if sensitivities are required
needs_sensitivities = False

# Logging
log_to_screen = True
log_filename = '/pstore/home/barraly/Codes/Gamma_0 impact on fitting/OHara initial conditions.csv'
log_csv = True
message_interval = 20
message_warm_up = 3

# Parallelisation
n_workers = pints.ParallelEvaluator.cpu_count()

#
# Stopping criteria
#

# Maximum iterations
max_iterations = 100000

# Maximum unchanged iterations
max_unchanged_iterations = 100
min_significant_change = 1

# Threshold value
threshold = None

# Post-run statistics
evaluations = None
iterations = None
time = None

# Iterations and function evaluations
iteration = 0
evaluations = 0

# Unchanged iterations count
unchanged_iterations = 0

# Choose method to evaluate
f = function
if needs_sensitivities:
    f = f.evaluateS1

# Create evaluator object
# For population based optimisers, don't use more workers than
# particles!
n_workers = min(n_workers, optimiser.population_size())
evaluator = pints.ParallelEvaluator(f, n_workers=n_workers)

# Keep track of best position and score
fbest = float('inf')
xbest = x0

# Internally we always minimise! Keep a 2nd value to show the user
fbest_user = fbest if minimising else -fbest

# Set up progress reporting
next_message = 0

# Start logging
if log_to_screen:
    # Show method
    print('Using ' + str(optimiser.name()))

    # Show parallelisation
    print('Running in parallel with ' + str(n_workers) +
          ' worker processes.')

# Show population size
pop_size = optimiser.population_size()
if log_to_screen:
    print('Population size: ' + str(pop_size))

# Set up logger
logger = pints.Logger()
logger.set_filename(log_filename, csv=log_csv)

# Add fields to log
max_iter_guess = max(max_iterations or 0, 10000)
max_eval_guess = max_iter_guess * pop_size
logger.add_counter('Iter.', max_value=max_iter_guess)
logger.add_counter('Eval.', max_value=max_eval_guess)
logger.add_float('Best')

#Log the values of parameters returning the best score
for param in range(sabs_pkpd.constants.n):
    logger.add_float(parameters_to_fit[param])

# Initialise logger
optimiser._log_init(logger)
logger.add_time('Time m:s')

# Start searching
timer = pints.Timer()

# Log first point
logger.log(0, 0, fbest_user, *xbest)
optimiser._log_write(logger)
logger.log(timer.time())

running = True
try:
    while running:
        # Get points
        xs = optimiser.ask()

        # Calculate scores
        fs = evaluator.evaluate(xs)

        # Perform iteration
        optimiser.tell(fs)

        # Check if new best found
        fnew = optimiser.fbest()
        if fnew < fbest:
            # Check if this counts as a significant change
            if np.abs(fnew - fbest) < min_significant_change:
                unchanged_iterations += 1
            else:
                unchanged_iterations = 0

            # Update best
            fbest = fnew
            xbest = optimiser.xbest()

            # Update user value of fbest
            fbest_user = fbest if minimising else -fbest
        else:
            unchanged_iterations += 1

        # Update evaluation count
        evaluations += len(fs)

        # Log state
        logger.log(iteration, evaluations, fbest_user, *xbest)
        optimiser._log_write(logger)
        logger.log(timer.time())

        # Update iteration count
        iteration += 1

        #
        # Check stopping criteria
        #

        # Maximum number of iterations
        if (max_iterations is not None and
                iteration >= max_iterations):
            running = False
            halt_message = ('Halting: Maximum number of iterations ('
                            + str(iteration) + ') reached.')

        # Maximum number of iterations without significant change
        halt = (max_unchanged_iterations is not None and
                unchanged_iterations >= max_unchanged_iterations)
        if halt:
            running = False
            halt_message = ('Halting: No significant change for ' +
                            str(unchanged_iterations) + ' iterations.')

        # Error in optimiser
        error = optimiser.stop()
        if error:   # pragma: no cover
            running = False
            halt_message = ('Halting: ' + str(error))
            

except (Exception, SystemExit, KeyboardInterrupt):  # pragma: no cover
    # Unexpected end!
    # Show last result and exit
    print('\n' + '-' * 40)
    print('Unexpected termination.')
    print('Current best score: ' + str(fbest))
    print('Current best position:')

    xbest = optimiser.xbest()

    for p in xbest:
        print(pints.strfloat(p))
    print('-' * 40)
    raise

# Stop timer
time = timer.time()

# Log final values and show halt message
logger.log(iteration, evaluations, fbest_user)
optimiser._log_write(logger)
logger.log(time)
print(halt_message)

# Save post-run statistics
evaluations = evaluations
iterations = iteration

# Inverse transform search parameters
xbest = optimiser.xbest()

print('/n /n Best point:/n')
print(xbest)

# Write the starting point to the last line in the format: 0, 0, 0, starting_point[0], ..., starting_point[-1], 0
with open(log_filename, 'a+', newline = '') as csvfile:
    mywriter = csv.writer(csvfile, delimiter=',')
    mywriter.writerow([0, 0, 0] + list(starting_point[:]) + [0])
    
    