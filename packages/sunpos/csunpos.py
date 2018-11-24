# encoding: utf-8
# module sunpos.csunpos
# from C:\Users\qinli\AppData\Local\Programs\Python\Python36\lib\site-packages\sunpos\csunpos.cp36-win_amd64.pyd
# by generator 1.145
# no doc

# imports
import builtins as __builtins__ # <module 'builtins' (built-in)>
import numpy as np # C:\Users\qinli\AppData\Local\Programs\Python\Python36\lib\site-packages\numpy\__init__.py
import datetime as __datetime


# Variables with simple values

__path__ = None

# functions

def arcdist(p0, p1, radians=False): # real signature unknown; restored from __doc__
    """
    arcdist(p0,p1,radians=False)
        Angular distance between azimuth,zenith pairs
        
        Parameters
        ----------
        p0 : array_like, shape (..., 2)
        p1 : array_like, shape (..., 2)
            p[...,0] = azimuth angles, p[...,1] = zenith angles
        radians : boolean (default False)
            If False, angles are in degrees, otherwise in radians
    
        Returns
        -------
        ad :  array_like, shape is broadcast(p0,p1).shape
            Arcdistances between corresponding pairs in p0,p1
            In degrees by default, in radians if radians=True
    """
    pass

def bibtex(*args, **kwargs): # real signature unknown
    pass

def citation(*args, **kwargs): # real signature unknown
    pass

def julian_day(dt): # real signature unknown; restored from __doc__
    """
    julian_day(dt)
        Convert UTC datetimes or UTC timestamps to Julian days
    
        Parameters
        ----------
        dt : array_like
            UTC datetime objects or UTC timestamps (as per datetime.utcfromtimestamp)
    
        Returns
        -------
        jd : ndarray
            datetimes converted to fractional Julian days
    """
    pass

def sunpos(dt, latitude, longitude, elevation, temperature=None, pressure=None, delta_t=0, radians=False): # real signature unknown; restored from __doc__
    """
    sunpos(dt, latitude, longitude, elevation, temperature=None, pressure=None, delta_t=0, radians=False)
        Compute the observed and topocentric coordinates of the sun as viewed at the given time and location.
    
        Parameters
        ----------
        dt : array_like
            UTC datetime objects or UTC timestamps (as per datetime.utcfromtimestamp) representing the times of observations
        latitude, longitude : array_like
            decimal degrees, positive for north of the equator and east of Greenwich
        elevation : array_like
            meters, relative to the WGS-84 ellipsoid
        temperature : array_like or None, optional
            celcius, default is 14.6 (global average in 2013)
        pressure : array_like or None, optional
            millibar, default is 1013 (global average in ??)
        delta_t : array_like, optional
            seconds, default is 0, difference between the earth's rotation time (TT) and universal time (UT)
        radians : {True, False}, optional
            return results in radians if True, degrees if False (default)
    
        Returns
        -------
        coords : ndarray, (...,5)
            The shape of the array is parameters broadcast together, plus a final dimension for the coordinates.
            coords[...,0] = observed azimuth angle, measured eastward from north
            coords[...,1] = observed zenith angle, measured down from vertical
            coords[...,2] = topocentric right ascension
            coords[...,3] = topocentric declination
            coords[...,4] = topocentric hour angle
    """
    pass

def sunpos_adh(dt, latitude, longitude, elevation, temperature=None, pressure=None, delta_t=0, radians=False): # real signature unknown; restored from __doc__
    """
    sunpos_adh(dt, latitude, longitude, elevation, temperature=None, pressure=None, delta_t=0, radians=False)
        Compute the right ascension, declination, and hour angles of the sun as viewed at the given time and location.
    
        Parameters
        ----------
        dt : array_like
            UTC datetime objects or UTC timestamps (as per datetime.utcfromtimestamp) representing the times of observations
        latitude, longitude : array_like
            decimal degrees, positive for north of the equator and east of Greenwich
        elevation : array_like
            meters, relative to the WGS-84 ellipsoid
        temperature : array_like or None, optional
            celcius, default is 14.6 (global average in 2013)
        pressure : array_like or None, optional
            millibar, default is 1013 (global average in ??)
        delta_t : array_like, optional
            seconds, default is 0, difference between the earth's rotation time (TT) and universal time (UT)
        radians : {True, False}, optional
            return results in radians if True, degrees if False (default)
    
        Returns
        -------
        coords : ndarray, (...,3)
            The shape of the array is parameters broadcast together, plus a final dimension for the coordinates.
            coords[...,0] = topocentric right ascension
            coords[...,1] = topocentric declination
            coords[...,2] = topocentric hour angle
    """
    pass

def sunpos_az(dt, latitude, longitude, elevation, temperature=None, pressure=None, delta_t=0, radians=False): # real signature unknown; restored from __doc__
    """
    sunpos_az(dt,latitude,longitude,elevation,temperature=None,pressure=None,delta_t=0,radians=False)
        Compute the azimuth and zenith angles of the sun as viewed at the given time and location.
    
        Parameters
        ----------
        dt : array_like
            UTC datetime objects or UTC timestamps (as per datetime.utcfromtimestamp) representing the times of observations
        latitude, longitude : array_like
            decimal degrees, positive for north of the equator and east of Greenwich
        elevation : array_like
            meters, relative to the WGS-84 ellipsoid
        temperature : array_like or None, optional
            celcius, default is 14.6 (global average in 2013)
        pressure : array_like or None, optional
            millibar, default is 1013 (global average in ??)
        delta_t : array_like, optional
            seconds, default is 0, difference between the earth's rotation time (TT) and universal time (UT)
        radians : {True, False}, optional
            return results in radians if True, degrees if False (default)
    
        Returns
        -------
        coords : ndarray, (...,2)
            The shape of the array is parameters broadcast together, plus a final dimension for the coordinates.
            coords[...,0] = observed azimuth angle, measured eastward from north
            coords[...,1] = observed zenith angle, measured down from vertical
    """
    pass

def __pyx_unpickle_Enum(*args, **kwargs): # real signature unknown
    pass

# classes

class datetime(__datetime.date):
    """
    datetime(year, month, day[, hour[, minute[, second[, microsecond[,tzinfo]]]]])
    
    The year, month and day arguments are required. tzinfo may be None, or an
    instance of a tzinfo subclass. The remaining arguments may be ints.
    """
    def astimezone(self, *args, **kwargs): # real signature unknown
        """ tz -> convert to local time in new timezone tz """
        pass

    @classmethod
    def combine(cls, *args, **kwargs): # real signature unknown
        """ date, time -> datetime with same date and time fields """
        pass

    def ctime(self): # real signature unknown; restored from __doc__
        """ Return ctime() style string. """
        pass

    def date(self, *args, **kwargs): # real signature unknown
        """ Return date object with same year, month and day. """
        pass

    def dst(self): # real signature unknown; restored from __doc__
        """ Return self.tzinfo.dst(self). """
        pass

    @classmethod
    def fromtimestamp(cls, *args, **kwargs): # real signature unknown
        """ timestamp[, tz] -> tz's local time from POSIX timestamp. """
        pass

    def isoformat(self, *args, **kwargs): # real signature unknown
        """
        [sep] -> string in ISO 8601 format, YYYY-MM-DDT[HH[:MM[:SS[.mmm[uuu]]]]][+HH:MM].
        sep is used to separate the year from the time, and defaults to 'T'.
        timespec specifies what components of the time to include (allowed values are 'auto', 'hours', 'minutes', 'seconds', 'milliseconds', and 'microseconds').
        """
        pass

    @classmethod
    def now(cls, *args, **kwargs): # real signature unknown
        """
        Returns new datetime object representing current time local to tz.
        
          tz
            Timezone object.
        
        If no tz is specified, uses local timezone.
        """
        pass

    def replace(self, *args, **kwargs): # real signature unknown
        """ Return datetime with new specified fields. """
        pass

    @classmethod
    def strptime(cls): # real signature unknown; restored from __doc__
        """ string, format -> new datetime parsed from a string (like time.strptime()). """
        pass

    def time(self, *args, **kwargs): # real signature unknown
        """ Return time object with same time but with tzinfo=None. """
        pass

    def timestamp(self, *args, **kwargs): # real signature unknown
        """ Return POSIX timestamp as float. """
        pass

    def timetuple(self, *args, **kwargs): # real signature unknown
        """ Return time tuple, compatible with time.localtime(). """
        pass

    def timetz(self, *args, **kwargs): # real signature unknown
        """ Return time object with same time and tzinfo. """
        pass

    def tzname(self): # real signature unknown; restored from __doc__
        """ Return self.tzinfo.tzname(self). """
        pass

    @classmethod
    def utcfromtimestamp(cls, *args, **kwargs): # real signature unknown
        """ Construct a naive UTC datetime from a POSIX timestamp. """
        pass

    @classmethod
    def utcnow(cls, *args, **kwargs): # real signature unknown
        """ Return a new datetime representing UTC day and time. """
        pass

    def utcoffset(self): # real signature unknown; restored from __doc__
        """ Return self.tzinfo.utcoffset(self). """
        pass

    def utctimetuple(self, *args, **kwargs): # real signature unknown
        """ Return UTC time tuple, compatible with time.localtime(). """
        pass

    def __add__(self, *args, **kwargs): # real signature unknown
        """ Return self+value. """
        pass

    def __eq__(self, *args, **kwargs): # real signature unknown
        """ Return self==value. """
        pass

    def __getattribute__(self, *args, **kwargs): # real signature unknown
        """ Return getattr(self, name). """
        pass

    def __ge__(self, *args, **kwargs): # real signature unknown
        """ Return self>=value. """
        pass

    def __gt__(self, *args, **kwargs): # real signature unknown
        """ Return self>value. """
        pass

    def __hash__(self, *args, **kwargs): # real signature unknown
        """ Return hash(self). """
        pass

    def __init__(self, year, month, day, hour=None, minute=None, second=None, microsecond=None, tzinfo=None): # real signature unknown; restored from __doc__
        pass

    def __le__(self, *args, **kwargs): # real signature unknown
        """ Return self<=value. """
        pass

    def __lt__(self, *args, **kwargs): # real signature unknown
        """ Return self<value. """
        pass

    @staticmethod # known case of __new__
    def __new__(*args, **kwargs): # real signature unknown
        """ Create and return a new object.  See help(type) for accurate signature. """
        pass

    def __ne__(self, *args, **kwargs): # real signature unknown
        """ Return self!=value. """
        pass

    def __radd__(self, *args, **kwargs): # real signature unknown
        """ Return value+self. """
        pass

    def __reduce_ex__(self, proto): # real signature unknown; restored from __doc__
        """ __reduce_ex__(proto) -> (cls, state) """
        pass

    def __reduce__(self): # real signature unknown; restored from __doc__
        """ __reduce__() -> (cls, state) """
        pass

    def __repr__(self, *args, **kwargs): # real signature unknown
        """ Return repr(self). """
        pass

    def __rsub__(self, *args, **kwargs): # real signature unknown
        """ Return value-self. """
        pass

    def __str__(self, *args, **kwargs): # real signature unknown
        """ Return str(self). """
        pass

    def __sub__(self, *args, **kwargs): # real signature unknown
        """ Return self-value. """
        pass

    fold = property(lambda self: object(), lambda self, v: None, lambda self: None)  # default

    hour = property(lambda self: object(), lambda self, v: None, lambda self: None)  # default

    microsecond = property(lambda self: object(), lambda self, v: None, lambda self: None)  # default

    minute = property(lambda self: object(), lambda self, v: None, lambda self: None)  # default

    second = property(lambda self: object(), lambda self, v: None, lambda self: None)  # default

    tzinfo = property(lambda self: object(), lambda self, v: None, lambda self: None)  # default


    max = None # (!) real value is ''
    min = None # (!) real value is ''
    resolution = None # (!) real value is ''


# variables with complex values

__loader__ = None # (!) real value is ''

__spec__ = None # (!) real value is ''

__test__ = {}

