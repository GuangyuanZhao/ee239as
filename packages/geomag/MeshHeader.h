#ifndef MESHHEADER_H
#define MESHHEADER_H

/*Mesh Functions*/
/*EMM specific constants, structures, and functions*/
typedef struct {
    double version; /* Model version */
    double epoch; /* Model year */

    int nalt; /* number of altitude layers */
    double *alt; /* array[nalt] */

    int *nlat; /* array[nalt] of number of rows in latitude */
    double *latres; /* array[nalt] containing the latitudinal resolution at each altitude */

    int **nlon; /* array[nalt][nlat] containing number of cells in each row */
    double **lonres; /* array[nalt][nlat] containing the longitudinal resolution in each row */

    double ****comp; /* variable array[nalt][nlat][nlon] containing data for x, y, and z vector components and derivatives */

} EMM_tmesh;


#define MIN_POLE_DIST   (0.001 * M_PI/180.0) /* Sperical coordinates not defined closer to poles */
#define USE_DERIVATIVES 1       /* Turn off only for testing */
#define INTERPOLATE_CUBIC 1     /* Turn off only for testing */
#define MESH_REAL_TYPE float    /* use 'float' to save space in binary file */

#define MESH_MAXLAT 1440            /* To check consistency of nlat parameter read from mesh file */

/* gridint return error codes */

#define EXIT_MESH_MEM_ALLOC_ERROR       1  /* Memory allocation failed */

#define EXIT_MESH_TEXT_FILE_NOT_FOUND   11 /* Input mesh file (text format) could not be opened */
#define EXIT_MESH_BIN_FILE_WRITE_ERROR  12 /* Output mesh file (binary format) could not be opened */
#define EXIT_MESH_BIN_FILE_NOT_FOUND    13 /* Binary mesh file could not be opened for reading */ 
#define EXIT_MESH_FILE_FORMAT_ERROR     14 /* Input mesh file format error. Last integer in binary file should have value 4711 */ 

#define EXIT_MESH_FILE_EPOCH_ERROR      20 /* Model year out of range */
#define EXIT_MESH_FILE_NALT_ERROR       21 /* Number of altitude layers out of range */
#define EXIT_MESH_FILE_NLAT_ERROR       22 /* Number of rows at this altitude is out of range */
#define EXIT_MESH_FILE_ALT_ERROR        23 /* Altitude out of range */
#define EXIT_MESH_FILE_NLON_ERROR       24 /* Number cells at this latitude and altitude is out of range */

#define WGS84_A 6378.1370       /* in km */
#define WGS84_B 6356.752314     /* in km */
#define WGS84_E sqrt(WGS84_A*WGS84_A - WGS84_B*WGS84_B)
#define WGS84_E2 (WGS84_E  * WGS84_E)
#define WGS84_E4 (WGS84_E2 * WGS84_E2)

extern void EMM_geodetic2geocentric(double phi, double h, double *latrad, double *r); /* phi is geodetic latitude in radian*/
extern void EMM_geocentric2geodetic_vec(double theta, double delta, double bx, double by, double bz,
        double *x, double *y, double *z);
extern void EMM_cart2sphere(double x, double y, double z, double *phi_rad, double *lat_rad, double *r);
extern void EMM_sphere2cart(double phi_rad, double lat_rad, double r, double *x, double *y, double *z);
extern void EMM_sphere2cart_vec(double phi_rad, double lat_rad, double vlat, double vphi, double vminusr,
        double *vx, double *vy, double *vz);
extern void EMM_cart2sphere_vec(double phi_rad, double lat_rad, double vx, double vy, double vz,
        double *vlat, double *vphi, double *vminusr);

int EMM_mesh_convert(int verbose, char infname[], char outfname[]);
int EMM_mesh_read(int verbose, char meshfname[], EMM_tmesh *mesh);
int EMM_Grid(MAGtype_CoordGeodetic minimum, MAGtype_CoordGeodetic maximum, double cord_step_size, double altitude_step_size, double time_step, MAGtype_Geoid *Geoid, MAGtype_Ellipsoid Ellip, MAGtype_Date StartDate, MAGtype_Date EndDate, int ElementOption, int PrintOption, char *OutputFile, EMM_tmesh mesh, EMM_tmesh meshSV);
int EMM_mesh_interpolate(int verbose, EMM_tmesh mesh, double lon, double lat, double alt,
        double geoc_lat, double *geoc_Bx, double *geoc_By, double *geoc_Bz);
int EMM_PointCalcFromMesh(MAGtype_CoordGeodetic CordGeo, MAGtype_CoordSpherical CordSph, MAGtype_Date UserDate, MAGtype_MagneticResults *MagResults, EMM_tmesh mesh, EMM_tmesh mesh_SV);

/*End of EMM specific constants, structures, and functions*/

#endif /*MESHHEADER_H*/
