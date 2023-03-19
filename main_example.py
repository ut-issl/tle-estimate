import BatchEstimator
from sgp4.api import Satrec
from sgp4.conveniences import jday_datetime
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
# Set times and dates

# GPS ON PASS no 190
# filename = "data/AOBC.AOBC_HK_COMPO.PRIVATE.NOT_REALTIME_#187ON_#190DR.csv"
# tle_celes = "1 55072U 23001BR  23040.15743944  .00016529  00000+0  91057-3 0  9998\n2 55072  97.5013 101.7243 0016943  94.7500 265.5667 15.13962877  5687"

# GPS ON PASS no 194
# filename = "data/AOBC.AOBC_HK_COMPO.PRIVATE.NOT_REALTIME_#194ON_#204DR.csv"
# tle_celes =  "1 55072U 23001BR  23053.93441113  .00010900  00000+0  59575-3 0  9997\n2 55072  97.4989 115.2981 0016335  49.7572 121.6134 15.14352866  7765"

# GPS ON PASS no 204
# filename = "data/AOBC.AOBC_HK_COMPO.PRIVATE.NOT_REALTIME_#204ON.csv"
# tle_celes =  "1 55072U 23001BR  23043.06547401  .00016320  00000+0  89684-3 0  9993\n2 55072  97.5006 104.5891 0016738  85.4079 274.9064 15.14056635  6121"

# GPS ON PASS no 220
# filename = "data/AOBC.AOBC_HK_COMPO.PRIVATE.NOT_REALTIME_#220ON.csv"
# tle_celes =  "1 55072U 23001BR  23047.55935536  .00018963  00000+0  10368-2 0  9998\n2 55072  97.5007 109.0179 0016616  70.8044 289.4984 15.14205762  6800"

# GPS ON PASS no 236
# filename = "data/AOBC.AOBC_HK_COMPO.PRIVATE.NOT_REALTIME_#236ON.csv"
# tle_celes =  "1 55072U 23001BR  23051.54355523  .00010900  00000+0  59657-3 0  9992\n2 55072  97.4978 112.9427 0016536  57.0774  48.6280 15.14300861  7400"

# GPS ON PASS no 261
# filename = "data/AOBC.AOBC_HK_COMPO.PRIVATE.NOT_REALTIME_#261ON.csv"
# tle_celes =  "1 55072U 23001BR  23058.47254928  .00012063  00000+0  65689-3 0  9992\n2 55072  97.4971 119.7701 0015750  35.7550  21.3429 15.14466278  8459"

# GPS ON PASS no 265
# filename = "data/AOBC.AOBC_HK_COMPO.PRIVATE.NOT_REALTIME_#265ON.csv"
# tle_celes =  "1 55072U 23001BR  23059.41090469  .00017025  00000+0  92399-3 0  9996\n2 55072  97.4970 120.6947 0015496  32.8711  97.0701 15.14508086  8597"

# GPS ON PASS no 274 - thrust maneouvre on
filename = "data/AOBC.AOBC_HK_COMPO.PRIVATE.NOT_REALTIME_#274ON.csv"
tle_celes =  "1 55072U 23001BR  23061.66111434  .00014480  00000+0  78517-3 0  9990\n2 55072  97.4979 122.9125 0015193  25.0534 126.1123 15.14581885  8938"

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
# code = batch_estimator.estimate_batch_orbit('ignore','j2')
if code < 0: print('Failed integration due to error (%d)' % code); exit()

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
pos_diff = np.zeros(batch_estimator.obs_num); vel_diff = np.zeros(num_states)
time = np.zeros(batch_estimator.obs_num)
for i in range(0,num_states):
    jd,fr = jday_datetime(batch_estimator.time[i])
    _,r_celes,v_celes = sat_celes.sgp4(jd, fr)
    _,r,v = sat.sgp4(jd, fr)
    pos_diff[i] = np.linalg.norm(np.array(r)-np.array(r_celes))
    vel_diff[i] = np.linalg.norm(np.array(v)-np.array(v_celes))
    time[i] = pd.Timedelta(batch_estimator.time[i] - batch_estimator.time[0]).total_seconds()
    
# Plot difference
plt.plot(time/60,pos_diff)
plt.ylabel('Position Difference [km]')
plt.xlabel('Time since epoch [min]')
plt.show()
plt.plot(time/60,vel_diff*1e3)
plt.ylabel('Velocity Difference [m/s]')
plt.xlabel('Time since epoch [min]')
plt.show()

