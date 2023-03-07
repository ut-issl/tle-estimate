"""@brief Python code that reads GPS dataset and calculates TLE """

##
# @mainpage Batch estimate from GPS to calculate TLE
#
# Author: Joshua Critchley-Marrows (ISSL/USYD) 2023-03-01
# Modified: Joshua Critchley-Marrows (ISSL/USYD) 2023-03-01

# Imports
import numpy as np
import pandas as pd
from numpy import float64
import pyatmos
import spiceypy as spice
from pyatmos.jb2008.spaceweather import read_sw_jb2008
import numdifftools as nd
import scipy as sp
import pyshtools as pysh
from sgp4.api import Satrec
from sgp4.conveniences import jday_datetime

# Global Constants
DEBUG = 1 # mode of operation (0 = normal, 1 = debug)
MAX_OBS_NUM = 10000 # maximum number of GPS observations
R_EARTH = 6378137 # Earth radius [m]
GRAV = 6.6743e-11 # gravitational constant [m3/kg/s2]
M_EARTH = 5.972e24 # mass of Earth [kg]
F_EARTH = 0.00335279499764807 # Earth flattening constant
MIN_VIS_SATS = 5 # minimum number of visible GPS satellites to include measurment
MODEL_GRAV = 'J2' 
MODEL_ATMOS = 'ignore'

class BatchEstimator:

    def __init__(self):
        """!@brief Initializes the program."""
        
        # Load the SPICE kernels
        spice.furnsh('batch_estimator.tm')
        
        # Update celestrack file
        self.swfile = pyatmos.download_sw_jb2008()
        
        # Assume BSTAR
        self.bstar = 0.59575e-3 # TLE btar [R_EARTH-1]
        
        # Read spherical harmonic file
        self.clm = pysh.datasets.Earth.EGM2008(2)
        
        if DEBUG:
            print("Initializing program.")
            
    def read_data_header(self,line):
        """!@brief Reads data header and returns to the list of parameters
        
        [time,GPS_TIME_week,GPS_TIME_msec,TLM_RECEIVED_TI,VISIBLE_SAT,
         POS_ECEF_X_m,POS_ECEF_Y_m,POS_ECEF_Z_m,VEL_ECEF_X_M/s,
         VEL_ECEF_Y_m/s,VEL_ECEF_Z_m/s]
         
         @param line    line containing the GPS data header
         
         @return header locations
        
        """
        header = line.split(',')
        header_loc = np.empty(11,dtype=int)
        header_loc.fill(-1)
        
        for i in range(len(header)):
            
            if header[i] == 'time':
                header_loc[0] = i
            elif header[i] == 'GPS_R.GPS_TIME_week':
                header_loc[1] = i
            elif header[i] == 'GPS_R.GPS_TIME_msec':
                header_loc[2] = i
            elif header[i] == 'GPS_R.TLM_RECEIVED_TI':
                header_loc[3] = i
            elif header[i] == 'GPS_R.VISIBLE_SAT':
                header_loc[4] = i
            elif header[i] == 'GPS_R.POS_ECEF_X_m':
                header_loc[5] = i
            elif header[i] == 'GPS_R.POS_ECEF_Y_m':
                header_loc[6] = i
            elif header[i] == 'GPS_R.POS_ECEF_Z_m':
                header_loc[7] = i
            elif header[i] == 'GPS_R.VEL_ECEF_X_m/s':
                header_loc[8] = i
            elif header[i] == 'GPS_R.VEL_ECEF_Y_m/s':
                header_loc[9] = i
            elif header[i] == 'GPS_R.VEL_ECEF_Z_m/s':
                header_loc[10] = i
        
        return header_loc
    
    def read_gps_data(self,filename):
        """!@brief Reads the GPS data set and retrieves GPS-related data.
                   Assigns to arrays
        
        @param    filename    File path to data set
        
        @return   Read success (1) or failure (0)
        """
        
        # Set constants
        GPS_EPOCH = pd.Timestamp('1980-01-06 00:00:00.0000',unit='ns',tz='utc')
        
        # Open the gps data file
        gps_data = open(filename)
        data_set = gps_data.readline()
        
        # Check the header locations
        header_loc = self.read_data_header(data_set)
        while np.sum(header_loc) == -11:
            data_set = gps_data.readline()
            header_loc = self.read_data_header(data_set)
            
            if ~data_set:
                return 0
            
        # Assign GPS dataset arrays
        self.rcv_time = np.empty((MAX_OBS_NUM,),dtype=object)
        self.gps_week = np.empty(MAX_OBS_NUM,dtype=float64)
        self.gps_msec = np.empty(MAX_OBS_NUM,dtype=float64)
        self.tlm_received = np.empty(MAX_OBS_NUM,dtype=float64)
        self.nsats = np.empty(MAX_OBS_NUM,dtype=int)
        self.states = np.empty((6,MAX_OBS_NUM),dtype=float64)
        self.obs_num = 0
        
        # Read GPS datasets from each line
        data_set = gps_data.readline()
        while data_set:
            
            data_set = data_set.split(',')
            
            # Time of observation            
            self.rcv_time[self.obs_num] = data_set[header_loc[0]]
                    
            # GPS Time of Week
            self.gps_week[self.obs_num] = data_set[header_loc[1]]
                    
            # GPS Milliseconds
            self.gps_msec[self.obs_num] = data_set[header_loc[2]]
                    
            # GPS Time Telemetry Receiver
            self.tlm_received[self.obs_num] = data_set[header_loc[3]]     
                            
            # GPS Number of Visible Satellites
            self.nsats[self.obs_num] = data_set[header_loc[4]]
                            
            # GPS Positions in ECEF [m]
            self.states[0,self.obs_num] = data_set[header_loc[5]]
            self.states[1,self.obs_num] = data_set[header_loc[6]]
            self.states[2,self.obs_num] = data_set[header_loc[7]]
                            
            # GPS Velocities in ECEF [m/s]
            self.states[3,self.obs_num] = data_set[header_loc[8]]
            self.states[4,self.obs_num] = data_set[header_loc[9]]
            self.states[5,self.obs_num] = data_set[header_loc[10]]
            
            # Move to the next line
            if self.nsats[self.obs_num] > MIN_VIS_SATS:
                self.obs_num = self.obs_num + 1
            data_set = gps_data.readline()
        
        # Reduce data sets to minimum number
        self.rcv_time = self.rcv_time[:self.obs_num]
        self.rcv_time = pd.to_datetime(self.rcv_time,yearfirst=True,unit='ns',utc=True)
        self.gps_week = self.gps_week[:self.obs_num]
        self.gps_msec = self.gps_msec[:self.obs_num]
        self.tlm_received = self.tlm_received[:self.obs_num]
        self.nsats = self.nsats[:self.obs_num]
        self.states = self.states[:,:self.obs_num]
        
        # Calculate gps weeks since epoch
        timeSinceGpsEpoch = (self.rcv_time[0] - GPS_EPOCH)
        weekPeriod = np.floor(timeSinceGpsEpoch.days / 7 / 2048)
        if weekPeriod*2048 > self.gps_week[0]:
            self.gps_week = self.gps_week + weekPeriod*2048
        
        # Convert GPS week and msec to measurement time
        self.time = GPS_EPOCH + pd.to_timedelta(self.gps_week,unit='W') + \
                    pd.to_timedelta(self.gps_msec,unit='ms')
        
        return 1
    
    def earth_gravity(self,time,r,model):
        """!@brief Determines the Earth's gravitational force and Jacobian matrix
        
        @para    time    Time value to derive the desired values (datetime)
        @para    r       Position of the  satellite in ECI (x,y,z) [m]
        @para    model   Treat Earth using spherical harmonics ('spherical') or as a point mass
        
        @return    a    Earth gravitation acceleration (a_x,a_y,a_z) [m/s2]
        @return    dadr Earth gravitation Jacobian in terms of position ((3,3))
        """
        
        J2 = 0.001082 # J2 perturbation non-dimensional
        J2 = J2*R_EARTH**2*GRAV*M_EARTH
        
        # Calculate distance
        rr = np.linalg.norm(r)
        
        if 'spherical' == model:
            
            # Transform to fixed frame
            et = spice.datetime2et(time)
            rot = spice.pxform('J2000','ITRF93',et)
            r_f = rot@r
            long,lat,_ = spice.recgeo(r_f,R_EARTH,F_EARTH)
            
            # Calculate gravitational acceleration
            a_sph = pysh.gravmag.MakeGravGridPoint(self.clm.coeffs,self.clm.gm,self.clm.r0,\
                                                   rr,lat*180/np.pi,long*180/np.pi)
            a = np.array((a_sph[0]*np.sin(a_sph[1])*np.cos(a_sph[2]),
                          a_sph[0]*np.sin(a_sph[1])*np.sin(a_sph[2]),
                          a_sph[0]*np.cos(a_sph[1])))
            dcm = [[-np.sin(lat)*np.cos(long), -np.sin(long), -np.cos(lat)*np.cos(long)],
                   [-np.sin(lat)*np.sin(long),  np.cos(long), -np.cos(lat)*np.sin(long)],
                   [np.cos(lat),                        0, -np.sin(lat)]]
            a = -rot.T@dcm@a

            # Calculate Jacobian
            dadr = self.spherical_gravity_jacobian(r_f)
            
        elif 'J2' == model:
            
            # Calculate terms
            x = r[0]; y = r[1]; z = r[2]
            term0 = 6*z**2 - 3/2*(x**2 + y**2)
            term1 = 3*z**2 - 9/2*(x**2 + y**2)

            # Calculate acceleration
            a = -GRAV*M_EARTH*r/rr**3
            a[0] = a[0] + J2*x/rr**7*term0
            a[1] = a[1] + J2*y/rr**7*term0
            a[2] = a[2] + J2*z/rr**7*term1
            
            # Calculate Jacobian
            dadr = -GRAV*M_EARTH*(np.ones((3,3))/rr**3 -\
                                   3*r.reshape((3,1))@r.reshape((3,1)).T/rr**5)
            dadr[0,0] = dadr[0,0] + J2*(rr**2-7*x**2)/rr**9*term0 - 3*J2*x**2/rr**7
            dadr[0,1] = dadr[0,1] - 7*J2*x*y/rr**9*term0 - 3*J2*x*y/rr**7
            dadr[0,2] = dadr[0,2] - 7*J2*x*z/rr**9*term0 + 12*J2*x*z/rr**7
            dadr[1,0] = dadr[1,0] - 7*J2*x*y/rr**9*term0 - 3*J2*x*y/rr**7
            dadr[1,1] = dadr[1,1] + J2*(rr**2-7*y**2)/rr**9*term0 - 3*J2*y**2/rr**7
            dadr[1,2] = dadr[1,2] - 7*J2*y*z/rr**9*term0 + 12*J2*y*z/rr**7
            dadr[2,0] = dadr[2,0] - 7*J2*x*z/rr**9*term1 - 9*J2*x*z/rr**7
            dadr[2,1] = dadr[2,1] - 7*J2*y*z/rr**9*term1 - 9*J2*y*z/rr**7
            dadr[2,2] = dadr[2,2] + J2*(rr**2-7*z**2)/rr**9*term1 + 6*J2*z**2/rr**7
            
            return a,dadr
        else:
        
            # Calculate acceleration
            a = -GRAV*M_EARTH*r/rr**3
            
            # Calculate Jacobian
            dadr = -GRAV*M_EARTH*(np.ones((3,3))/rr**3 -\
                                   3*r.reshape((3,1))@r.reshape((3,1)).T/rr**5)
        
        return a,dadr
    
    def find_atmospheric_density(self,time,r,model):
        """!@brief Finds the atmospheric density at a certain position
        
        @para    time    Time value to derive the desired values (datetime)
        @para    r       Position of the  satellite in ECI (x,y,z) [m]
        @para    model   Use either JB2008 ('jb2008') or exponential model
        
        @return    rho   Atmospheric density [kg/m3]
        """
                
        # Calculate transformation matrix
        et = spice.datetime2et(time)
        rot = spice.pxform('J2000','ITRF93',et)
        r = rot@r
        
        # Calculate lat, long, alt
        long,lat,alt = spice.recgeo(r.reshape(-1),R_EARTH,F_EARTH)
        lat = lat*180/np.pi; long = long*180/np.pi; alt = alt/1000
        
        # Estimate atmospheric density (assume spherical Earth) (use exponential model)
        if model == 'jb2008':
            swdata = read_sw_jb2008(self.swfile)
            jb08 = pyatmos.jb2008(time,(lat,long,alt),swdata)
            rho = jb08.rho
        else:
            expom = pyatmos.expo(alt)
            rho = expom.rho

        return rho        
    
    def atmospheric_drag(self,time,r,v,model):
        """!@brief Determines the Earth's atmospheric drag and Jacobians
        
        @para    time    Time value to derive the desired values (datetime)
        @para    r       Position of the satellite in ECI (x,y,z) [m]
        @para    v       Velocity of the  atellite in ECI (v_x,v_y,v_z) [m]
        @para    model   Treat Earth using spherical harmonics ('spherical') or as a point mass
        
        @return    a    Earth drag acceleration (a_x,a_y,a_z) [m/s2]
        @return    dadr Earth drag Jacobian in terms of position ((3,3))
        @return    dadv Earth drag Jacobian in terms of velocity ((3,3))
        """
        
        # Define constants
        RHO_0 = 0.1570 # reference air density [kg/m2/R_EARTH]
        
        if model == 'ignore':
            return np.zeros(3),np.zeros((3,3)), np.zeros((3,3))
        
        # Calculate gravitational density
        rho = self.find_atmospheric_density(time,r.reshape((3,1)),model)
        
        # Relative velocity of spacecraft relative to Earth atmosphere (assume relative to
        # Earth rotation)
        et = spice.datetime2et(time)
        rot = spice.sxform('ITRF93','J2000',et)
        _,omega = spice.xf2rav(rot)
        v_r = v - np.cross(omega,r)
        
        # Calculate atmospheric drag
        a = -rho/(RHO_0*R_EARTH)*self.bstar*v_r*np.linalg.norm(v_r)
        
        # Estimate drag derivative
        fun = lambda rr: self.find_atmospheric_density(time,rr,'expo')
        drhodr = nd.Jacobian(fun)
        drhodr = drhodr(r.reshape((3,1)))
        
        # Calculate Jacobian
        X = [[0.,-omega[2],omega[1]],
             [omega[2],0.,-omega[0]],
             [-omega[1],omega[0],0.]]
        v_r = v_r.reshape((3,1))
        dadv = -rho/(RHO_0*R_EARTH)*(v_r@v_r.T/np.linalg.norm(v_r) + \
                                     np.linalg.norm(v_r)*np.ones((3,1)))
        dadr = -rho/(RHO_0*R_EARTH)*np.linalg.norm(v_r)*v_r@drhodr - \
                    dadv@X
        
        return a,dadr,dadv
    
    def estimate_batch_orbit(self):
        """!@brief  Calculates the batch orbit estimates using SVD from a set of position and
                    velocity measurements
                    Reference:
                    Gill, E., Montenbruck, O. (2000). Satellite orbits: models, 
                    methods, and applications. Germany: Springer Berlin Heidelberg.
            
            @return success result or return error code of integration tool
        """
        
        # Assume initial position and velocity as GPS estimate
        et = spice.datetime2et(self.time[0])
        rot = spice.sxform('ITRF93','J2000',et)
        x0_original = rot@self.states[:,0]
        
        # Set up estimation matrices
        # self.obs_num = 1
        H = np.empty((6*self.obs_num,6))
        dz = np.empty((6*self.obs_num,1))
        self.x0 = x0_original
        dx = np.ones((6,1))

        # Build the estimation operation
        iter_count = 0
        while np.linalg.norm(dx) > 1e-5 and iter_count < 10:
            for i in range(self.obs_num):
                
                # Compute state and transformation
                et = spice.datetime2et(self.time[i])
                rot = spice.sxform('ITRF93','J2000',et)
    
                # Retrieve measurement
                z = rot@np.reshape(self.states[:,i],(6,1))
                
                # Set current time
                dt = pd.Timedelta(self.time[i] - self.time[0]).total_seconds()
                
                # Calculate propagated state x approximate by numerical integration
                def xfun(t,x,t0):
                    time = t0 + pd.to_timedelta(t, unit='sec')
                    a_e,_ = self.earth_gravity(time,x[0:3],MODEL_GRAV)
                    a_atm,_,_ = self.atmospheric_drag(time,x[0:3],x[3:6],MODEL_ATMOS)
                    a = a_e + a_atm
                    return np.concatenate((x[3:6], a)) 
                integ = sp.integrate.ode(xfun).set_integrator('vode')
                integ.set_initial_value(self.x0,0).set_f_params(self.time[0])
                x = integ.integrate(dt).reshape((6,1))
                if integ.successful() != 1:   
                    return integ.get_return_code()
                
                # Calculate acting forces
                _,dadr_e = self.earth_gravity(self.time[0],x0_original[0:3],MODEL_GRAV)
                _,dadr_atm,dadv_atm = self.atmospheric_drag(self.time[i],x0_original[0:3],x0_original[3:6],MODEL_ATMOS)
                
                # Calculate F matrix
                F = np.block([[np.zeros((3,3)), np.eye(3)],
                              [dadr_e + dadr_atm, dadv_atm]])
                
                # Calculate Phi approximate by numerical integration
                def Phifun(t,Phi,F):
                    Phi = Phi.reshape((6,6))
                    dPhi = F@Phi
                    return dPhi.reshape(-1)
                dt = pd.Timedelta(self.time[i] - self.time[0]).total_seconds()
                Phi0 = np.eye(6).reshape(-1)
                integ = sp.integrate.ode(Phifun).set_integrator('vode')
                integ.set_initial_value(Phi0,0).set_f_params(F)
                Phi = integ.integrate(dt).reshape((6,6))
                if integ.successful() != 1:
                    return integ.get_return_code()
                
                # Contribute this to H and z vectors
                dz[i*6:(i*6+6),:] = z - x
                H[i*6:(i*6+6),:] = Phi
                
            
            # Calculate SVD estimation
            U,D,V = sp.linalg.svd(H,full_matrices=False)
            D = np.diag(D)
            dx = V.T@sp.linalg.inv(D)@U.T@dz
            self.x0 = self.x0 + dx.reshape(-1)
 
            iter_count = iter_count + 1

        
        self.t0 = self.time[0]
        
        return 1
    
    def estimate_batch_orbit_sgp4(self):
        """!@brief  Calculates the batch orbit estimates using SVD from a set of position and
                    velocity measurements, where the SGP4 propagator is used.
            
            @return success result or return error code of integration tool
        """
        
        # Set optimisation parameters
        tol = 1e-3
        delta_perc = 0.001
        max_iter = 100
        
        # Define solver function that will need to equal zero for Kepler 
        # to be exact to sgp4 TLE
        def sgpsolver(para,tle,time):
            
            # Set basic TLE
            s_tle = tle.split('\n')[0]
            t_tle = tle.split('\n')[1]
            t_tle = t_tle[:8] + '%8.4f' % para[2] + t_tle[16:]
            t_tle =  t_tle[:17] + '%8.4f' % para[3] + t_tle[25:]
            ecc = '%.7f' % para[1]
            t_tle = t_tle[:26] + '' + ecc.lstrip('0').lstrip('.') + t_tle[33:]
            t_tle = t_tle[:34] + '%8.4f' % para[4] + t_tle[42:]
            t_tle = t_tle[:43] + '%8.4f' % para[5] + t_tle[51:]
            t_tle = t_tle[:52] + '%11.8f' % para[0] + t_tle[63:]
            
            # Set up satellite object
            sat = Satrec.twoline2rv(s_tle, t_tle)
            
            # Run SGP4 propagation with time difference as zero
            jd,fr = jday_datetime(time)
            _,r,v = sat.sgp4(jd,fr)
            
            # Calculate tolerance 
            x = np.zeros(6)
            x[0:3] = np.array(r)*1e3; x[3:6] = np.array(v)*1e3
            return x
                
        # Initialise Kepler parameters
        et = spice.datetime2et(self.time[0])
        rot = spice.sxform('ITRF93','J2000',et)
        self.x0 = rot@self.states[:,0]
        self.t0 = self.time[0]
        self.state_to_kepler()
        self.write_to_tle()
        self.state_to_kepler_sgp4()
        self.write_to_tle()
    
        # Set up estimation matrices
        H = np.empty((6*self.obs_num,6))
        dz = np.empty((6*self.obs_num,1))

        # Build the estimation operation
        para = self.para.copy()
        deltay = 1.0
        iter_count = 0
        while np.linalg.norm(deltay) > tol and max_iter > iter_count:
            for i in range(self.obs_num):
                
                # Transform measured state to the ECI/J2000 frame
                et = spice.datetime2et(self.time[i])
                rot = spice.sxform('ITRF93','J2000',et)
                x = rot@np.reshape(self.states[:,i],(6,1))
                
                # Calculate propagated state x approximate by numerical integration
                M = np.zeros((6,6))
                f = sgpsolver(para,self.tle_str,self.time[i])
                for j in range(0,6): 
                    delta_para = para.copy()
                    delta_para[j] = para[j] + delta_perc*para[j]
                    deltaf = sgpsolver(delta_para,self.tle_str,self.time[i])
                    M[:,j] = (deltaf - f)/(delta_perc*para[j])
        
                # Contribute this to H and z vectors
                dz[i*6:(i*6+6),:] = x - f.reshape((6,1))
                H[i*6:(i*6+6),:] = M
                
            # Calculate SVD estimation
            U,D,V = sp.linalg.svd(H,full_matrices=False)
            D = np.diag(D)
            deltay = V.T@sp.linalg.inv(D)@U.T@dz
            para = para + deltay.reshape(-1)
            
            # Check if all terms are within bounds
            if para[1] > 1 or para[1] < 0: para[1] = self.para[1]
            para[2] = para[2] % 360; para[3] = para[3] % 360
            para[4] = para[4] % 360; para[5] = para[5] % 360
 
            iter_count = iter_count + 1

        
        self.t0 = self.time[0]
        self.para = para
        self.write_to_tle()
        
        return 1
    
    def state_to_kepler(self):
        """!@brief Transform ECI position and velocity to Kepler orbital parameters.
        """
        
        mu = GRAV*M_EARTH
        
        # Calculate position and velocity norm
        r = self.x0[0:3]
        v = self.x0[3:6]
        rr = np.linalg.norm(r)
        vv = np.linalg.norm(v)
        
        # Calculate angular momentum and node vector
        h = np.cross(r,v)
        k = [0.,0.,1.]
        N = np.cross(k,h)
        
        # Calculate eccentricity vector
        ee = ((vv**2 - mu/rr)*r - np.dot(r,v)*v)/mu
        ecc = np.linalg.norm(ee) 
        
        # Calculate energy
        energy = vv**2/2 - mu/rr
        
        # Determine semimajor axis and mean motion
        a_major = -mu/2/energy
        n_motion = np.sqrt(mu/a_major**3)
        n_motion = n_motion/2/np.pi*60*60*24
        
        # Determine inclination
        inc = np.arccos(h[2]/np.linalg.norm(h))
        inc = inc*180./np.pi
        
        # Determine right ascension of ascending node
        raan = np.arccos(N[0]/np.linalg.norm(N))
        if N[1] < 0: raan = 2*np.pi - raan
        raan = raan*180/np.pi
        
        # Determine argument of perigee
        aop = np.arccos(np.dot(N,ee)/np.linalg.norm(N)/ecc)
        if ee[2] < 0: aop = 2*np.pi - aop
        aop = aop*180./np.pi
        
        # Determine true anomaly
        true_anom = np.arccos(np.dot(ee,r)/ecc/rr)
        if np.dot(r,v) < 0: true_anom = 2*np.pi - true_anom
        
        # Determine eccentric anomaly
        E_anom = 2*np.arctan2(((1+ecc)/(1-ecc))**(-1/2)*np.sin(true_anom/2),\
                       np.cos(true_anom/2))
            
        # Determine mean anomaly at epoch
        M_anom = E_anom - ecc*np.sin(E_anom)
        M_anom = M_anom*180/np.pi
        
        # Set parameter
        # [n_motion, ecc, inc, raan, aop, M_anom]
        self.para = [n_motion, ecc, inc, raan, aop, M_anom]
          
        return
    
    def state_to_kepler_sgp4(self):
        """!@brief Transform ECI position and velocity to Kepler orbital parameters by the SGP4
                   definition.
                   Reference: Lee, Byoung-Sun & Park, Jae-Woo. (2003). Estimation of the SGP4 Drag Term
                   From Two Osculating Orbit States. Journal of Astronomy and Space Sciences. 20. 11-20. 
                   10.5140/JASS.2003.20.1.011. 
        """
        
        # Set constants
        tol = 1e-5
        delta_perc = 0.001
        max_iter = 100
        
        # Define solver function that will need to equal zero for Kepler 
        # to be exact to sgp4 TLE
        def sgpsolver(para,tle,time):
            
            # Set basic TLE
            s_tle = tle.split('\n')[0]
            t_tle = tle.split('\n')[1]
            t_tle = t_tle[:8] + '%8.4f' % para[2] + t_tle[16:]
            t_tle =  t_tle[:17] + '%8.4f' % para[3] + t_tle[25:]
            ecc = '%.7f' % para[1]
            t_tle = t_tle[:26] + '' + ecc.lstrip('0').lstrip('.') + t_tle[33:]
            t_tle = t_tle[:34] + '%8.4f' % para[4] + t_tle[42:]
            t_tle = t_tle[:43] + '%8.4f' % para[5] + t_tle[51:]
            t_tle = t_tle[:52] + '%11.8f' % para[0] + t_tle[63:]
            
            # Set up satellite object
            jd,fr = jday_datetime(time)
            sat = Satrec.twoline2rv(s_tle, t_tle)
            
            # Run SGP4 propagation with time difference as zero
            _,r,v = sat.sgp4(jd,fr)
            
            # Calculate tolerance 
            x = np.zeros(6)
            x[0:3] = np.array(r)*1e3; x[3:6] = np.array(v)*1e3
            return x
        
        # Set up iterative solver loop
        para = self.para.copy()
        deltay = 1.0
        iter_count = 0
        while np.linalg.norm(deltay) > tol and max_iter > iter_count:
            
            # Build Jacobian using forward difference quotient
            M = np.zeros((6,6))
            f = sgpsolver(para,self.tle_str,self.t0)
            for i in range(0,6): 
                delta_para = para.copy()
                delta_para[i] = para[i] + delta_perc*para[i]
                deltaf = sgpsolver(delta_para,self.tle_str,self.t0)
                M[:,i] = (deltaf - f)/(delta_perc*para[i])
        
            # Calculate iteration
            Minv = np.linalg.inv(M)
            deltax = self.x0 - f 
            deltay = Minv@deltax
            para = para + deltay
            
            # Check if all terms are within bounds
            if para[1] > 1 or para[1] < 0: para[1] = self.para[1]
            para[2] = para[2] % 360; para[3] = para[3] % 360
            para[4] = para[4] % 360; para[5] = para[5] % 360
            
            iter_count = iter_count + 1
        
        # Set parameter
        # [n_motion, ecc, inc, raan, aop, M_anom]
        self.para = para
          
        return
    
    
    def read_tle(self,tle_str):
        """!@brief Read and stores the parameters extracted from a TLE string.
        
            @para    tle    TLE string
        
            @return    Sucessful operation boolean
        """
        
        # Split tle lines
        tle = tle_str.split('\n')
        
        # Check tle is first row
        if tle[0][0:1] != '1': return 0
        
        # Read first line contents
        self.tle_sat_cat_num = np.int(tle[0][2:7])
        self.tle_classification = tle[0][7:8]
        self.tle_int_desig = tle[0][9:17]
        year = pd.to_datetime(tle[0][18:20],format="%y")
        doy = pd.to_timedelta(np.double(tle[0][20:32]),unit='D')
        self.t0 = year + doy - pd.Timedelta(days=1)
        self.tle_ndot = np.double(tle[0][33:43])
        self.tle_nddot = np.double('.' + tle[0][44:50].strip() \
                                   + 'e' + tle[0][50:52])
        self.bstar = np.double(tle[0][53:54] + '.' + tle[0][54:59].strip() \
                                   + 'e' + tle[0][59:61])
        self.tle_eph_type = tle[0][62:63]
        self.ele_set_no = np.int(tle[0][64:68])
                        
        # Check tle is second row
        if tle[1][0:1] != '2': return 0
        
        # Read second line contents
        inc = np.double(tle[1][8:16])
        raan = np.double(tle[1][17:25])
        ecc = np.double('.' + tle[1][26:33].strip())
        aop = np.double(tle[1][33:42])
        M_anom = np.double(tle[1][43:51])
        n_motion = np.double(tle[1][52:63])
        self.tle_rev_num = np.int(tle[1][63:68])
        self.para = [n_motion, ecc, inc, raan, aop, M_anom]
        
        # Store string parameter in class
        self.tle_str = tle_str
        
        return 1
    
    def convert_sci(self,n):
        """!@brief Convert scientific notation-based number to TLE format 'XXXXX...X+X'
        
        @para    n    Number of significant figures
        
        @return    String contain TLE format based scientific notation-based number
        """
        sign = '-' if self.tle_nddot < 0 else ' '
        if n == 0.0: return ' 00000+0'
        ppow = np.floor(np.log10(n))
        if ppow < -9: return ' 00000+0'
        pow_str = '%1d' % (np.abs(ppow)-1)
        pow_sign = '-' if ppow < 0 else ' '
        dec = (np.abs(n)/10**ppow)*10**-1
        dec_str = '%5.5f' % dec
        
        return sign + dec_str.lstrip('0').lstrip('.') + pow_sign + pow_str
    
    def calc_checksum(self,line):
        """!@brief Calculate checksum of the line by the summation of all digits, modulo 10
        
        @para    line    Line to calculate checksum over
        
        @return   checksum
        """
        
        checksum = 0
        for ch in line:
            if ch.isdigit() == True: checksum += np.int(ch)
            elif ch == '-': checksum += 1
        
        return checksum % 10
    
    def write_to_tle(self):
        """!@brief Write stored parameters to a tle string
        
            @return    TLE string
        """
        
        # Write first line of tle
        tle = '1'
        tle = tle + ' %05d' % self.tle_sat_cat_num
        tle = tle + self.tle_classification
        tle = tle + ' ' + self.tle_int_desig
        tle = tle + ' ' + self.t0.strftime('%y')
        doy = self.t0.dayofyear + self.t0.hour / 24.0 + self.t0.minute / (24.0*60.0) + \
                (self.t0.second + np.double(self.t0.strftime('.%f'))) / (24.0*60.0*60.0) 
        tle = tle + '%012.8f' % doy
        sign = ' -' if self.tle_ndot < 0 else '  '
        ndot = '%8.8f' % self.tle_ndot
        tle = tle + sign + ndot.lstrip('0')
        tle = tle + ' ' + self.convert_sci(self.tle_nddot)
        tle = tle + ' ' + self.convert_sci(self.bstar)
        tle = tle + ' ' + self.tle_eph_type
        tle = tle + ' %4d' % self.ele_set_no
        checksum = self.calc_checksum(tle)
        tle = tle + '%1d' % checksum
        tle = tle +'\n'
        
        # Write second line of tle
        tle_2 = '2'
        tle_2 = tle_2 + ' %05d' % self.tle_sat_cat_num
        tle_2 = tle_2 + ' %8.4f' % self.para[2]
        tle_2 = tle_2 + ' %8.4f' % self.para[3]
        ecc = '%.7f' % self.para[1]
        tle_2 = tle_2 + ' ' + ecc.lstrip('0').lstrip('.')
        tle_2 = tle_2 + ' %8.4f' % self.para[4]
        tle_2 = tle_2 + ' %8.4f' % self.para[5]
        tle_2 = tle_2 + ' %11.8f' % self.para[0]
        rev_num = '%5d' % self.tle_rev_num
        tle_2 = tle_2 + rev_num.ljust(5)
        checksum = self.calc_checksum(tle_2)
        tle_2 = tle_2 + '%1d' % checksum
        tle_2 = tle_2 +'\n'
        tle = tle + tle_2
        
        self.tle_str = tle
        
        return tle
        
    def spherical_gravity_jacobian(self,x):
        """!@brief Calculate the spherical harmonic Jacobian
            Reference: Lundberg, J.B., & Schutz, B.E. (1988). Recursion formulas of Legendre functions 
                       for use with nonsingular geopotential models. Journal of Guidance Control and 
                       Dynamics, 11, 31-38.
        
        @para    x    Position of the satellite in ECEF coordinates (x,y,z) [m]
        
        @return       Spherical harmonic Jacobian matrix
        """
        
        # Read coefficients
        Cnm = self.clm.coeffs[0,:,:]
        Snm = self.clm.coeffs[1,:,:]
        if self.clm.normalization == '4pi':
            nor = 4*np.pi
        else:
            nor = 1
        degree = self.clm.lmax
        Cnm = Cnm*nor; Snm = Snm*nor
        
        # Derive new spherical coordinates
        r = np.linalg.norm(x);
        s = x[0]/r; t = x[1]/r; u = x[2]/r
    
        # Calculate r_m and i_m
        # [r_0 r_1 ... r_n], [i_0 i_1 ... i_n]
        r_m = np.ones(degree+1); i_m = np.zeros(degree+1)
        r_m[1] = s; i_m[1] = t
        for m in range(2,degree+1):
            r_m[m] = s*r_m[m-1] - t*i_m[m-1]
            i_m[m] = s*i_m[m-1] + t*r_m[m-1]
    
        # Calculate Anm terms
        # [[A_00  A_01  ... A_0n2 ]
        #  [A_10  A_11  ... A_1n2 ]
        #  [...   ...   ... ...   ] 
        #  [A_n20 A_n21 ... A_n2n2]]
        Anm = np.zeros((degree+3,degree+3))
        cosphi = np.sqrt(1 - x[2]**2/r**2)
        Anm[1,1] = np.sqrt(3)*cosphi
        for n in range(2,degree+3):
            for m in range(1,n+1):
                if m == n:
                    Anm[n,m] = cosphi*np.sqrt((2*n+1)/2/n)*Anm[n-1,m-1]
                elif n == 2:
                    Anm[n,m] = u*np.sqrt((2*n+1)*(2*n-1)/(n-m)/(n+m))*Anm[n-1,m]
                else:
                    Anm[n,m] = u*np.sqrt((2*n+1)*(2*n-1)/(n-m)/(n+m))*Anm[n-1,m] - \
                        np.sqrt((2*n+1)*(n-m-1)*(n+m-1)/(2*n-3)/(n+m)/(n-m))*Anm[n-2,m]
    
        # Calculate rho terms
        # [rho0 rho1 rho2 ... rho_n+2]
        rho = np.zeros(degree+3);
        rho[0] = self.clm.gm/r;
        for n in range(1,degree+3):
            rho[n] = self.clm.r0/r*rho[n-1];
    
        # Calculate Dnm to Hnm terms
        # note: same shape as Anm
        Dnm = np.zeros((degree+1,degree+1)); Enm = np.zeros((degree+1,degree+1))
        Fnm = np.zeros((degree+1,degree+1)); Gnm = np.zeros((degree+1,degree+1))
        Hnm = np.zeros((degree+1,degree+1))
        for n in range(0,degree+1):
            for m in range(0,n+1):
                Nnm = np.sqrt(np.math.factorial(n-m)*np.math.factorial(2*n+1)/np.math.factorial(n+m));
                Dnm[n,m] = Cnm[n,m]/Nnm*r_m[m] + Snm[n,m]/Nnm*i_m[m];
                if m > 0:
                    Enm[n,m] = Cnm[n,m]/Nnm*r_m[m-1] + Snm[n,m]/Nnm*i_m[m-1];
                    Fnm[n,m] = Snm[n,m]/Nnm*r_m[m-1] - Cnm[n,m]/Nnm*i_m[m-1];
                if m > 1:
                    Gnm[n,m] = Cnm[n,m]/Nnm*r_m[m-2] + Snm[n,m]/Nnm*i_m[m-2];
                    Hnm[n,m] = Snm[n,m]/Nnm*r_m[m-2] - Cnm[n,m]/Nnm*i_m[m-2];
    
        # Calculate a terms
        a_11 = 0; a_12 = 0; a_13 = 0; a_14 = 0; a_23 = 0; a_24 = 0
        a_33 = 0; a_34 = 0; a_44 = 0; a_4 = 0
        for n  in range(0,degree+1):  
            for m in range(0,n+1):
    
                cnm1 = np.sqrt((n-m)*(n+m+1))
                cn1m1 = np.sqrt((n+m+2)*(n+m+1)/(2*n+3)/(2*n+2))
                cn1m2 = cn1m1*np.sqrt((n-m)*(n+m+3))
                cnm2 = cnm1*np.sqrt((n-m+1)*(n+m+2))
                cn2m2 = cn1m1*np.sqrt((n+m+4)*(n+m+3)/(2*n+5)/(2*n+4))
    
                factor = rho[n+2]/self.clm.r0**2
                a_11 = a_11 + factor*m*(m-1)*Anm[n,m]*Gnm[n,m]
                a_12 = a_12 + factor*m*(m-1)*Anm[n,m]*Hnm[n,m]
                a_13 = a_13 + factor*m*cnm1*Anm[n,m+1]*Enm[n,m]
                a_14 = a_14 - factor*m*cn1m1*Anm[n+1,m+1]*Enm[n,m]
                a_23 = a_23 + factor*m*cnm1*Anm[n,m+1]*Fnm[n,m]
                a_24 = a_24 - factor*m*cn1m1*Anm[n+1,m+1]*Fnm[n,m]
                a_33 = a_33 + factor*cnm2*Anm[n,m+2]*Dnm[n,m]
                a_34 = a_34 - factor*cn1m2*Anm[n+1,m+2]*Dnm[n,m]
                a_44 = a_44 + factor*cn2m2*Anm[n+2,m+2]*Dnm[n,m]
    
                a_4 = a_4 - (rho[n+1]/self.clm.r0)*cn1m1*Anm[n+1,m+1]*Dnm[n,m]
    
        # Calculate final terms
        C = np.zeros((3,3))
        C[0,0] = a_11 + s**2*a_44 + a_4/r + 2*s*a_14
        C[0,1] = a_12 + s*t*a_44 + s*a_24 + t*a_14; C[1,0] = C[0,1]
        C[0,2] = a_13 + s*u*a_44 + s*a_34 + u*a_14; C[2,0] = C[0,2]
        C[1,1] = -a_11 + t**2*a_44 + a_4/r + 2*t*a_24
        C[1,2] = a_23 + t*u*a_44 + u*a_24 + t*a_34; C[2,1] = C[1,2]
        C[2,2] = a_33 + u**2*a_44 + a_4/r + 2*u*a_34
    
        # Assign nonspherical gravity gradient
        return C
        
        
        
        