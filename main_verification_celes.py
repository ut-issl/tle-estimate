import BatchEstimator
from sgp4.api import Satrec
from sgp4.conveniences import jday_datetime
import datetime
import numpy as np
import matplotlib.pyplot as plt

# Set times and dates
filename = "data/AOBC.AOBC_HK_COMPO.PRIVATE.NOT_REALTIME_#187ON_#190DR.csv"
# filename = "data/AOBC.AOBC_HK_COMPO.PRIVATE.NOT_REALTIME_#194ON_#204DR.csv"
tle_celes = "1 55072U 23001BR  23040.15743944  .00016529  00000+0  91057-3 0  9998\n2 55072  97.5013 101.7243 0016943  94.7500 265.5667 15.13962877  5687\n"
# tle_celes = "1 55072U 23001BR  23041.67757084  .00016707  00000+0  91904-3 0  9995\n2 55072  97.5007 103.2219 0016831  89.9008 270.4152 15.14012539  5917\n"

date = datetime.datetime.strptime('2023-02-11 12:00:00','%Y-%m-%d %H:%M:%S')

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
code = batch_estimator.estimate_batch_orbit()
if code < 0: print('Failed integration due to error (%d)' % code); exit()

# Convert estimated initial state to Kepler parameters
batch_estimator.state_to_kepler()

# Write calculated parameters to TLE
tle = batch_estimator.write_to_tle()
print(tle_celes)
print(tle)

# Compare to sgp4
s = tle.split('\n')[0]
t = tle.split('\n')[1]
sat_celes = Satrec.twoline2rv(s_celes, t_celes)
sat = Satrec.twoline2rv(s, t)
jd,fr = jday_datetime(date)
e_celes,r_celes,v_celes = sat_celes.sgp4(jd, fr)
e,r,v = sat.sgp4(jd, fr)
print(e_celes,r_celes,v_celes)
print(e,r,v)
