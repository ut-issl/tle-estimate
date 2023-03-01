import BatchEstimator
from sgp4.api import Satrec
from sgp4.conveniences import jday_datetime
import datetime

# Set times and dates
filename = "data/AOBC.AOBC_HK_COMPO.PRIVATE.NOT_REALTIME_#187ON_#190DR.csv"
# tle_celes =  "1 55072U 23001BR  23053.93441113  .00010900  00000+0  59575-3 0  9997\n2 55072  97.4989 115.2981 0016335  49.7572 121.6134 15.14352866  7765"
tle_celes = "1 55072U 23001BR  23040.15743944  .00016529  00000+0  91057-3 0  9998\n2 55072  97.5013 101.7243 0016943  94.7500 265.5667 15.13962877  5687"
tle =        "1 55072U 23001BR  23040.05457176  .00010900  00000+0  59575-3 0  9991\n2 55072  97.5815 101.4117 0052089 144.8391  14.1792 15.04649224  7766"
date = datetime.datetime.strptime('2023-03-01 00:00:00','%Y-%d-%m %H:%M:%S')

# Read Celestrack based TLE code
s_celes = tle_celes.split('\n')[0]
t_celes = tle_celes.split('\n')[1]
print(tle_celes)

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

