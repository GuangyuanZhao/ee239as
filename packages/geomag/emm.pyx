
cdef extern from "EGM9615.h":
    float *GeoidHeightBuffer

cdef extern from "GeomagnetismHeader.h":
    #types    
    ctypedef struct MAGtype_MagneticModel:
        pass

    ctypedef struct MAGtype_Ellipsoid:
        pass

    ctypedef struct MAGtype_CoordGeodetic:
        double lon 'lambda'
        double lat 'phi'
        double HeightAboveEllipsoid
        double HeightAboveGeoid
        int UseGeoid

    ctypedef struct MAGtype_CoordSpherical:
        pass

    ctypedef struct MAGtype_Date:
        double DecimalYear

    ctypedef struct MAGtype_LegendreFunction:
        pass

    ctypedef struct MAGtype_MagneticResults:
        double Bx # /* North */
        double By # /* East */
        double Bz # /* Down */

    ctypedef struct MAGtype_SphericalHarmonicVariables:
        pass

    ctypedef struct MAGtype_GeoMagneticElements:
        double Decl # /* 1. Angle between the magnetic field vector and true north, positive east*/
        double Incl # /*2. Angle between the magnetic field vector and the horizontal plane, positive down*/
        double F # /*3. Magnetic Field Strength*/
        double H # /*4. Horizontal Magnetic Field Strength*/
        double X # /*5. Northern component of the magnetic field vector*/
        double Y # /*6. Eastern component of the magnetic field vector*/
        double Z # /*7. Downward component of the magnetic field vector*/
        double GV # /*8. The Grid Variation*/
        double Decldot # /*9. Yearly Rate of change in declination*/
        double Incldot # /*10. Yearly Rate of change in inclination*/
        double Fdot # /*11. Yearly rate of change in Magnetic field strength*/
        double Hdot # /*12. Yearly rate of change in horizontal field strength*/
        double Xdot # /*13. Yearly rate of change in the northern component*/
        double Ydot # /*14. Yearly rate of change in the eastern component*/
        double Zdot # /*15. Yearly rate of change in the downward component*/
        double GVdot #

    ctypedef struct MAGtype_Geoid:
        float *GeoidHeightBuffer #
        int Geoid_Initialized # /* indicates successful initialization */
        int UseGeoid # /*Is the Geoid being used?*/

    ctypedef struct MAGtype_Gradient:
        pass

    ctypedef struct MAGtype_CoordGeodeticStr:
        pass

    ctypedef struct MAGtype_UTMParameters:
        pass

    #functions
    bint MAG_robustReadMagneticModel_Large(char *filename, char *filenameSV, MAGtype_MagneticModel **MagneticModel)
    bint MAG_SetDefaults(MAGtype_Ellipsoid *Ellip, MAGtype_Geoid *Geoid)
    bint MAG_GeodeticToSpherical(MAGtype_Ellipsoid Ellip, MAGtype_CoordGeodetic CoordGeodetic, MAGtype_CoordSpherical *CoordSpherical)
    bint MAG_ConvertGeoidToEllipsoidHeight(MAGtype_CoordGeodetic *CoordGeodetic, MAGtype_Geoid *Geoid)
    bint MAG_CalculateGeoMagneticElements(MAGtype_MagneticResults *MagneticResultsGeo, MAGtype_GeoMagneticElements *GeoMagneticElements)
    bint MAG_CalculateSecularVariationElements(MAGtype_MagneticResults MagneticVariation, MAGtype_GeoMagneticElements *MagneticElements)
    void MAG_AssignMagneticModelCoeffs(MAGtype_MagneticModel *ass, MAGtype_MagneticModel *src, int nMax, int nMaxSecVar)
    MAGtype_MagneticModel *MAG_AllocateModelMemory(int numTerms)

cdef extern from "MeshHeader.h":
    #types
    ctypedef struct EMM_tmesh:
        pass
    #error codes
    cdef enum:
        EXIT_MESH_MEM_ALLOC_ERROR      = 1  # Memory allocation failed 

        EXIT_MESH_TEXT_FILE_NOT_FOUND  = 11 # Input mesh file (text format) could not be opened
        EXIT_MESH_BIN_FILE_WRITE_ERROR = 12 # Output mesh file (binary format) could not be opened
        EXIT_MESH_BIN_FILE_NOT_FOUND   = 13 # Binary mesh file could not be opened for reading 
        EXIT_MESH_FILE_FORMAT_ERROR    = 14 # Input mesh file format error. Last integer in binary file should have value 4711 

        EXIT_MESH_FILE_EPOCH_ERROR     = 20 # Model year out of range
        EXIT_MESH_FILE_NALT_ERROR      = 21 # Number of altitude layers out of range
        EXIT_MESH_FILE_NLAT_ERROR      = 22 # Number of rows at this altitude is out of range
        EXIT_MESH_FILE_ALT_ERROR       = 23 # Altitude out of range
        EXIT_MESH_FILE_NLON_ERROR      = 24 # Number cells at this latitude and altitude is out of range
    #functions
    int EMM_mesh_read(int verbose, char *meshfname, EMM_tmesh *mesh)
    int EMM_PointCalcFromMesh(MAGtype_CoordGeodetic CordGeo, MAGtype_CoordSpherical CordSph, MAGtype_Date UserDate, MAGtype_MagneticResults *MagResults, EMM_tmesh mesh, EMM_tmesh mesh_SV)

cdef class GeoMagneticElements:
    """GeoMagneticElements encapsulates the geomagnetic field paramters at a location.
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
    cdef public double Decl, Incl, F, H, X, Y, Z, GV, Decldot, Incldot, Fdot, Hdot, Xdot, Ydot, Zdot, GVdot
    def __cinit__(self):
        Decl = Incl = F = H = X = Y = Z = GV = Decldot = Incldot = Fdot = Hdot = Xdot = Ydot = Zdot = GVdot = 0

    cdef _setup(self, MAGtype_GeoMagneticElements elements):
        self.Decl = elements.Decl
        self.Incl = elements.Incl
        self.F = elements.F
        self.H = elements.H
        self.X = elements.X
        self.Y = elements.Y
        self.Z = elements.Z
        self.GV = elements.GV
        self.Decldot = elements.Decldot
        self.Incldot = elements.Incldot
        self.Fdot = elements.Fdot
        self.Hdot = elements.Hdot
        self.Xdot = elements.Xdot
        self.Ydot = elements.Ydot
        self.Zdot = elements.Zdot
        self.GVdot = elements.GVdot

from libc.math cimport cos, sin, sqrt, atan2, asin
from time import gmtime
from os import path

cdef double d2r(double d):
    return d*0.017453292519943295
cdef double r2d(double r):
    return r*57.29577951308232

cdef class EMMBase:
    cdef MAGtype_Geoid _geoid
    cdef MAGtype_Ellipsoid _ellip

    def __cinit__(self):
        #initialize the ellipsoid and geoid
        if not MAG_SetDefaults(&self._ellip, &self._geoid):
            raise RuntimeError('Error setting default ellipsoid and geoid parameters.')
        self._geoid.GeoidHeightBuffer = GeoidHeightBuffer
        self._geoid.Geoid_Initialized = 1

    cdef tuple norm_lat_lon(self, lat, lon):
        """Wrap out of bounds lat, lon into the valid range"""
        if lat < -90 or lat > 90:
            #convert to cartesian and back
            x = cos(d2r(lon))*cos(d2r(lat))
            y = sin(d2r(lon))*cos(d2r(lat))
            z = sin(d2r(lat))
            r = sqrt(x**2 + y**2 + z**2)
            lon = r2d(atan2(y,x)) % 360
            lat = r2d(asin(z/r))
        elif lon < 0 or lon > 360:
            lon = lon % 360
        return lat,lon

    @staticmethod
    def citation():
        return "Chulliat, A., P. Alken, M. Nair, A. Woods, and S. Maus, 2015, The Enhanced Magnetic Model 2015-2020, National Centers for Environmental Information, NOAA. doi: 10.7289/V56971HV."
        
    @staticmethod
    def bibtex():
        return '@article{RN121,\n    author = {Chulliat, A and Alken, P and Nair, M and Woods, A and Maus, S},\n    title = {The Enhanced Magnetic Model 2015-2020},\n    journal = {National Centers for Environmental Information, NOAA},\n    DOI = {10.7289/V56971HV},\n    year = {2015},\n    type = {Journal Article}\n}'


    def decimal_year(self, t):
        """decimal_year(t)
        convert a time.struct_time, datetime.date, datetime.datetime, or UTC timestamp into a decimal year"""
        if hasattr(t,'timetuple'): #date/datetime to struct_time:
            tt = t.timetuple()
        elif hasattr(t,'tm_year') and hasattr(t,'tm_yday'): #already a struct_time
            tt = t
        else: #seconds to struct_time
            tt = gmtime(t)
        return tt.tm_year + (tt.tm_yday-1)/365 #leap years?

    @staticmethod
    cdef _EMM_check(int err):
        """Check error codes"""
        if err == EXIT_MESH_MEM_ALLOC_ERROR:
            raise MemoryError()
        elif err == EXIT_MESH_TEXT_FILE_NOT_FOUND:
            raise RuntimeError('Text file not found.')
        elif err == EXIT_MESH_BIN_FILE_WRITE_ERROR:
            raise RuntimeError('Binary file write error.')
        elif err == EXIT_MESH_BIN_FILE_NOT_FOUND:
            raise RuntimeError('Binary file not found.')
        elif err == EXIT_MESH_FILE_FORMAT_ERROR:
            raise RuntimeError('File format error.')
        elif err == EXIT_MESH_FILE_EPOCH_ERROR:
            raise RuntimeError('File epoch error (Model year out of range)')
        elif err == EXIT_MESH_FILE_NALT_ERROR:
            raise RuntimeError('File error: number of altitude layers out of range.')
        elif err == EXIT_MESH_FILE_NLAT_ERROR:
            raise RuntimeError('File error: number of latitude rows out of range.')
        elif err == EXIT_MESH_FILE_ALT_ERROR:
            raise RuntimeError('File error: altitude out of range.')
        elif err == EXIT_MESH_FILE_NLON_ERROR:
            raise RuntimeError('File error: number of longitude cells at latitude and altitude out of range.')

##cdef class EMMSph(EMMBase):
##    """EMMSph(cof_dir, first_year, last_year)
##        cof_dir: path to directory containing coefficient files
##        first_year, last_year: inclusive year range for loading coefficient files
##
##    This class wraps NOAA's Enhanced Magnetic Model (EMM) spherical harmonic routines
##    These routines have a lower memory footprint than EMMMesh, but will use more CPU time.
##    """


cdef class EMMMesh(EMMBase):
    """EMMMesh(str mesh_fname, str svmesh_fname, bool delay_load=False)
        mesh_fname: filename of EMM static mesh
        secmesh_fname: filename of EMM secular variation mesh
        delay_load: if True, meshes will not be loaded until load() or a function requiring them is called
    This class wraps NOAA's Enhanced Magnetic Model (EMM) Mesh routines.
    These routines use less CPU time than EMMSph, but have a larger memory footprint.
    """
    cdef EMM_tmesh _mesh
    cdef EMM_tmesh _mesh_sv
    cdef bytes _mfname
    cdef bytes _smfname
    cdef bint _mloaded
    cdef bint _smloaded

    def __cinit__(self, str mesh_fname, str secmesh_fname, bint delay_load = False):
        self._mfname = bytes(mesh_fname,'UTF-8')
        self._smfname = bytes(secmesh_fname,'UTF-8')
        self._mloaded = 0
        self._smloaded = 0
        if not delay_load:
            self._load_c()

    cdef MAGtype_GeoMagneticElements _compute_field_c(self, MAGtype_CoordGeodetic pos, MAGtype_CoordSpherical cs, MAGtype_Date date, bint compute_change = False):
        #local variables:
        cdef MAGtype_MagneticResults magResults, magVariation
        cdef MAGtype_GeoMagneticElements elements
        #load meshes if necessary:
        self._load_c()
        #compute magnetic results
        EMM_PointCalcFromMesh(pos, cs, date, &magResults, self._mesh, self._mesh_sv)
        #convert to geo-magnetic elements
        MAG_CalculateGeoMagneticElements(&magResults, &elements)
        #compute secular variation
        if compute_change:
            date.DecimalYear += 1
            EMM_PointCalcFromMesh(pos, cs, date, &magVariation, self._mesh, self._mesh_sv)
            magVariation.Bx -= magResults.Bx
            magVariation.By -= magResults.By
            magVariation.Bz -= magResults.Bz
            MAG_CalculateSecularVariationElements(magVariation, &elements)
        return elements

    cdef MAGtype_GeoMagneticElements compute_field_c(self, double lat, double lon, double height, double year, bint geodetic = True, bint compute_change = False):
        #local variables
        cdef MAGtype_Date date
        cdef MAGtype_CoordGeodetic pos
        cdef MAGtype_CoordSpherical cs

        #put year into a date struct:
        date.DecimalYear = year

        #convert coordinates
        pos.lat = lat
        pos.lon = lon
        if geodetic:
            pos.HeightAboveGeoid = height
            pos.UseGeoid = 1
            self._geoid.UseGeoid = 1
            MAG_ConvertGeoidToEllipsoidHeight(&pos, &self._geoid)
        else:
            pos.UseGeoid = 0
            self._geoid.UseGeoid = 0
            pos.HeightAboveGeoid = height
            pos.HeightAboveEllipsoid = height

        MAG_GeodeticToSpherical(self._ellip, pos, &cs)

        return self._compute_field_c(pos, cs, date, compute_change)

    cdef void _load_c(self):
        if not self._mloaded:
            EMMMesh._EMM_check(EMM_mesh_read(0, self._mfname, &self._mesh))
            self._mloaded = 1
        if not self._smloaded:
            EMMMesh._EMM_check(EMM_mesh_read(0, self._smfname, &self._mesh_sv))
            self._smloaded = 1

    def load(self,mesh=None,secmesh=None):
        """load(mesh=None,secmesh=None)
            Load the specified mesh files, or if delay_load=True was specified in the constructor, load those files.
            Paramters:
                mesh, secmesh - filenames for mesh or secular variation mesh. Both are optional
            No return value.
        """
        if mesh is not None: 
            self._mfname = bytes(mesh,'UTF-8')
            self._mloaded = False
        if secmesh is not None:
            self._smfname = bytes(secmesh, 'UTF-8')
            self._smloaded = False
        self._load_c()

    def is_loaded(self):
        """is_loaded()
        Return True if mesh files have been loaded
        """
        return self._mloaded and self._smloaded
        
    def compute_field(self, lat, lon, height, year, geodetic = True, compute_change = False):
        """compute_field(lat, lon, height, year, geodetic = True, compute_change = False)
            Parameters:
                lat - latitude, in degrees
                lon - longitude, in degrees
                height - height above EGM96 mean sea level, or WGS-84 ellipsoid if geodetic = False
                year - date, in decimal years
                geodetic - if true, use EGM96 mean sea level as reference for height, otherwise use WGS-84 ellipsoid
                compute_change - if true, compute secular variation of magnetic field (rate of chage per year)
            Returns a GeoMagneticElements object with the results.
        """
        lat,lon = self.norm_lat_lon(lat,lon)
        elements = GeoMagneticElements()
        elements._setup(self.compute_field_c(lat,lon,height,year,geodetic,compute_change))
        return elements

    def declination(self, lat, lon, height, year, geodetic = True):
        """declination(lat, lon, height, year, geodetic = True)
        Angle (deg) between the magnetic field vector and true north, positive east."""
        return self.compute_field(lat,lon,height,year,geodetic).Decl

    def true_to_magnetic(self, head, lat, lon, height, year, geodetic = True):
        """true_to_magnetic(head, lat, lon, height, year, geodetic = True)
        Convert a heading (deg) from true north to magnetic north (subtract declination)"""
        return head - self.declination(lat,lon,height,year,geodetic)

    def magnetic_to_true(self,head,lat,lon,height,year,geodetic=True):
        """magnetic_to_true(head, lat, lon, height, year, geodetic = True)
        Convert a heading (deg) from magnetic north to true north (add declination)"""
        return head + self.declination(lat,lon,height,year,geodetic)
    