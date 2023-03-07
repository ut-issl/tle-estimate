import BatchEstimator
from sgp4.api import Satrec
from sgp4.conveniences import jday_datetime
import datetime
import matplotlib.pyplot as plt
import numpy as np

# Set times and dates
filename = "data/AOBC.AOBC_HK_COMPO.PRIVATE.NOT_REALTIME_#187ON_#190DR.csv"
# filename = "data/AOBC.AOBC_HK_COMPO.PRIVATE.NOT_REALTIME_#194ON_#204DR.csv"
# tle_celes =  "1 55072U 23001BR  23053.93441113  .00010900  00000+0  59575-3 0  9997\n2 55072  97.4989 115.2981 0016335  49.7572 121.6134 15.14352866  7765"
tle_celes = "1 55072U 23001BR  23040.15743944  .00016529  00000+0  91057-3 0  9998\n2 55072  97.5013 101.7243 0016943  94.7500 265.5667 15.13962877  5687"
date = datetime.datetime.strptime('2023-02-11 00:00:00','%Y-%m-%d %H:%M:%S')

# Read Celestrack based TLE code
s_celes = tle_celes.split('\n')[0]
t_celes = tle_celes.split('\n')[1]

# Initialise batch estimator
batch_estimator = BatchEstimator.BatchEstimator()

# Read estimated TLE line
code = batch_estimator.read_tle(tle_celes)
if code == 0: print('Incorrect tle string format.'); exit()

# Read GPS data file
code = batch_estimator.read_gps_data(filename)
if code == 0: print('Unsuccessful in reading GPS data containing file'); exit()

# Estimate orbit based on GPS data
code = batch_estimator.estimate_batch_orbit_sgp4()
if code < 0: print('Failed integration due to error (%d)' % code); exit()

# Convert estimated initial state to Kepler parameters
batch_estimator.state_to_kepler()
# batch_estimator.state_to_kepler_sgp4()

# Write calculated parameters to TLE
tle = batch_estimator.write_to_tle()
print(tle_celes)
print(tle)

# Compare using sgp4
s = tle.split('\n')[0]
t = tle.split('\n')[1]
sat_celes = Satrec.twoline2rv(s_celes, t_celes)
sat = Satrec.twoline2rv(s, t)
num_states = batch_estimator.obs_num
pos_error = np.zeros(batch_estimator.obs_num); vel_error = np.zeros(num_states)
for i in range(0,num_states):
    jd,fr = jday_datetime(batch_estimator.time[i])
    _,r_celes,v_celes = sat_celes.sgp4(jd, fr)
    _,r,v = sat.sgp4(jd, fr)
    pos_error[i] = np.linalg.norm(np.array(r)-np.array(r_celes))
    vel_error[i] = np.linalg.norm(np.array(v)-np.array(v_celes))
    
# Plot error
plt.plot(pos_error)
plt.ylabel('Error [km]')
plt.xlabel('Observation Number')
plt.show()

