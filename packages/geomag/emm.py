# encoding: utf-8
# module geomag.emm
# from C:\Users\qinli\AppData\Local\Programs\Python\Python36\lib\site-packages\geomag\emm.cp36-win_amd64.pyd
# by generator 1.145
# no doc

# imports
import builtins as __builtins__ # <module 'builtins' (built-in)>
import ntpath as path # C:\Users\qinli\AppData\Local\Programs\Python\Python36\lib\ntpath.py
from time import gmtime


# Variables with simple values

__path__ = None

# no functions
# classes

class EMMBase(object):
    # no doc
    def bibtex(self, *args, **kwargs): # real signature unknown
        pass

    def citation(self, *args, **kwargs): # real signature unknown
        pass

    def decimal_year(self, t): # real signature unknown; restored from __doc__
        """
        decimal_year(t)
                convert a time.struct_time, datetime.date, datetime.datetime, or UTC timestamp into a decimal year
        """
        pass

    def __init__(self, *args, **kwargs): # real signature unknown
        pass

    @staticmethod # known case of __new__
    def __new__(*args, **kwargs): # real signature unknown
        """ Create and return a new object.  See help(type) for accurate signature. """
        pass

    def __reduce__(self, *args, **kwargs): # real signature unknown
        pass

    def __setstate__(self, *args, **kwargs): # real signature unknown
        pass

    __pyx_vtable__ = None # (!) real value is ''


class EMMMesh(EMMBase):
    """
    EMMMesh(str mesh_fname, str svmesh_fname, bool delay_load=False)
            mesh_fname: filename of EMM static mesh
            secmesh_fname: filename of EMM secular variation mesh
            delay_load: if True, meshes will not be loaded until load() or a function requiring them is called
        This class wraps NOAA's Enhanced Magnetic Model (EMM) Mesh routines.
        These routines use less CPU time than EMMSph, but have a larger memory footprint.
    """
    def compute_field(self, lat, lon, height, year, geodetic=True, compute_change=False): # real signature unknown; restored from __doc__
        """
        compute_field(lat, lon, height, year, geodetic = True, compute_change = False)
                    Parameters:
                        lat - latitude, in degrees
                        lon - longitude, in degrees
                        height - height above EGM96 mean sea level, or WGS-84 ellipsoid if geodetic = False
                        year - date, in decimal years
                        geodetic - if true, use EGM96 mean sea level as reference for height, otherwise use WGS-84 ellipsoid
                        compute_change - if true, compute secular variation of magnetic field (rate of chage per year)
                    Returns a GeoMagneticElements object with the results.
        """
        pass

    def declination(self, lat, lon, height, year, geodetic=True): # real signature unknown; restored from __doc__
        """
        declination(lat, lon, height, year, geodetic = True)
                Angle (deg) between the magnetic field vector and true north, positive east.
        """
        pass

    def is_loaded(self): # real signature unknown; restored from __doc__
        """
        is_loaded()
                Return True if mesh files have been loaded
        """
        pass

    def load(self, mesh=None, secmesh=None): # real signature unknown; restored from __doc__
        """
        load(mesh=None,secmesh=None)
                    Load the specified mesh files, or if delay_load=True was specified in the constructor, load those files.
                    Paramters:
                        mesh, secmesh - filenames for mesh or secular variation mesh. Both are optional
                    No return value.
        """
        pass

    def magnetic_to_true(self, head, lat, lon, height, year, geodetic=True): # real signature unknown; restored from __doc__
        """
        magnetic_to_true(head, lat, lon, height, year, geodetic = True)
                Convert a heading (deg) from magnetic north to true north (add declination)
        """
        pass

    def true_to_magnetic(self, head, lat, lon, height, year, geodetic=True): # real signature unknown; restored from __doc__
        """
        true_to_magnetic(head, lat, lon, height, year, geodetic = True)
                Convert a heading (deg) from true north to magnetic north (subtract declination)
        """
        pass

    def __init__(self, str_mesh_fname, str_svmesh_fname, bool_delay_load=False): # real signature unknown; restored from __doc__
        pass

    @staticmethod # known case of __new__
    def __new__(*args, **kwargs): # real signature unknown
        """ Create and return a new object.  See help(type) for accurate signature. """
        pass

    def __reduce__(self, *args, **kwargs): # real signature unknown
        pass

    def __setstate__(self, *args, **kwargs): # real signature unknown
        pass

    __pyx_vtable__ = None # (!) real value is ''


class GeoMagneticElements(object):
    """
    GeoMagneticElements encapsulates the geomagnetic field paramters at a location.
        members: 
          Decl = Angle between the magnetic field vector and true north, positive east.
          Incl = Angle between the magnetic field vector and the horizontal plane, positive down.
          F = Magnetic field strength.
          H = Horizontal magnetic field strength.
          X = Northern component of the magnetic field vector.
          Y = Eastern component of the magnetic field vector.
          Z = Downward component of the magnetic field vector.
          GV = The grid variation.
          Decldot, Incldot, Fdot, Hdot, Xdot, Ydot, Zdot, and GVdot: the change per year of each of the above.
    """
    def __init__(self, *args, **kwargs): # real signature unknown
        pass

    @staticmethod # known case of __new__
    def __new__(*args, **kwargs): # real signature unknown
        """ Create and return a new object.  See help(type) for accurate signature. """
        pass

    def __reduce__(self, *args, **kwargs): # real signature unknown
        pass

    def __setstate__(self, *args, **kwargs): # real signature unknown
        pass

    Decl = property(lambda self: object(), lambda self, v: None, lambda self: None)  # default

    Decldot = property(lambda self: object(), lambda self, v: None, lambda self: None)  # default

    F = property(lambda self: object(), lambda self, v: None, lambda self: None)  # default

    Fdot = property(lambda self: object(), lambda self, v: None, lambda self: None)  # default

    GV = property(lambda self: object(), lambda self, v: None, lambda self: None)  # default

    GVdot = property(lambda self: object(), lambda self, v: None, lambda self: None)  # default

    H = property(lambda self: object(), lambda self, v: None, lambda self: None)  # default

    Hdot = property(lambda self: object(), lambda self, v: None, lambda self: None)  # default

    Incl = property(lambda self: object(), lambda self, v: None, lambda self: None)  # default

    Incldot = property(lambda self: object(), lambda self, v: None, lambda self: None)  # default

    X = property(lambda self: object(), lambda self, v: None, lambda self: None)  # default

    Xdot = property(lambda self: object(), lambda self, v: None, lambda self: None)  # default

    Y = property(lambda self: object(), lambda self, v: None, lambda self: None)  # default

    Ydot = property(lambda self: object(), lambda self, v: None, lambda self: None)  # default

    Z = property(lambda self: object(), lambda self, v: None, lambda self: None)  # default

    Zdot = property(lambda self: object(), lambda self, v: None, lambda self: None)  # default


    __pyx_vtable__ = None # (!) real value is ''


# variables with complex values

__loader__ = None # (!) real value is ''

__spec__ = None # (!) real value is ''

__test__ = {}

