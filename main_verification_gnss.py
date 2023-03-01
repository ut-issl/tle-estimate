import BatchEstimator
from sgp4.api import Satrec
from sgp4.conveniences import jday_datetime
import datetime
import numpy as np
import matplotlib.pyplot as plt

# Set times and dates
filename = "data/AOBC.AOBC_HK_COMPO.PRIVATE.NOT_REALTIME_#187ON_#190DR.csv"
# tle_celes =  "1 55072U 23001BR  23053.93441113  .00010900  00000+0  59575-3 0  9997\n2 55072  97.4989 115.2981 0016335  49.7572 121.6134 15.14352866  7765"
tle_celes = "1 55072U 23001BR  23040.15743944  .00016529  00000+0  91057-3 0  9998\n2 55072  97.5013 101.7243 0016943  94.7500 265.5667 15.13962877  5687"
date = datetime.datetime.strptime('2023-03-01 00:00:00','%Y-%d-%m %H:%M:%S')

########### Compare GNSS measurements to batch estimate ##############

# Initialise batch estimator
batch_estimator = BatchEstimator.BatchEstimator()

# Read estimated TLE line
code = batch_estimator.read_tle(tle_celes)
if code == 0: print('Incorrect tle string format.'); exit()

# Read GPS data file
code = batch_estimator.read_gps_data(filename)
if code == 0: print('Unsuccessful in reading GPS data containing file'); exit()

# Estimate orbit based on GPS data
code = batch_estimator.estimate_batch_orbit()
if code < 0: print('Failed integration due to error (%d)' % code); exit()

# Convert estimated initial state to Kepler parameters
batch_estimator.state_to_kepler()

# Write calculated parameters to TLE
tle = batch_estimator.write_to_tle()

# Calculate orbit via SGP4
s = tle.split('\n')[0]
t = tle.split('\n')[1]
sat = Satrec.twoline2rv(s, t)
jd,fr = jday_datetime(date)
e_batch,r_batch,v_batch = sat.sgp4(jd, fr)

# Calculate TLE for individual measurements
states = batch_estimator.states
obs_num = batch_estimator.obs_num
times = batch_estimator.time
batch_estimator.obs_num = 1
pos_error = np.zeros(obs_num); vel_error = np.zeros(obs_num)
for i in range(0,obs_num):
    
    batch_estimator.states = states[:,i].reshape((6,1))
    batch_estimator.time = times[i:i+1]
    
    # Estimate orbit based on GPS data
    code = batch_estimator.estimate_batch_orbit()
    if code < 0: print('Failed integration due to error (%d)' % code); exit()
    
    # Convert estimated initial state to Kepler parameters
    batch_estimator.state_to_kepler()
    
    # Write calculated parameters to TLE adn derive sgp4
    tle = batch_estimator.write_to_tle()
    s = tle.split('\n')[0]
    t = tle.split('\n')[1]
    sat = Satrec.twoline2rv(s, t)
    
    # Calculate position and velocity errors
    e,r,v = sat.sgp4(jd, fr)
    pos_error[i] = np.linalg.norm(np.array(r_batch)-np.array(r))
    vel_error[i] = np.linalg.norm(np.array(v_batch)-np.array(v))

# Plot error
plt.plot(pos_error)
plt.ylabel('Error [km]')
plt.xlabel('Observation Number')
plt.show()

h = 2

