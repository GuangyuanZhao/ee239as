#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>

#include "GeomagnetismHeader.h"
#include "MeshHeader.h"

/*MESH STATIC DECLARATIONS*/

static int EMM_find_alt_index(int verbose, EMM_tmesh mesh, double alt, int *ialt0, int *ialt1,
        double *alt0, double *alt1, double *altfac0, double *altfac1);
static void EMM_find_lat_index(EMM_tmesh mesh, int ialt, double lat, int *ilat0, int *ilat1,
        double *lat0, double *lat1, double *latfac0, double *latfac1);
static void EMM_find_lon_index(EMM_tmesh mesh, int ialt, int ilat, double lon, int *ilon0, int *ilon1,
        double *lon0, double *lon1, double *lonfac0, double *lonfac1);
static double EMM_cell(EMM_tmesh mesh, int ialt, int ilat, int ilon, int icomp);
static double EMM_dlon_cell(EMM_tmesh mesh, int ialt, int ilat, int ilon, int icomp);
static double EMM_dlat_cell(EMM_tmesh mesh, int ialt, int ilat, int ilon, int icomp);
static double EMM_dalt_cell(EMM_tmesh mesh, int ialt, int ilat, int ilon, int icomp);
static double EMM_interpolate_quadratic(double x1, double y1, double dx_y1, double x2,
        double y2, double dx_y2, double x);
static double EMM_interpolate_cubic(double x1, double y1, double dx_y1, double x2,
        double y2, double dx_y2, double x);
static double EMM_interpolate(double x1, double y1, double dx_y1, double x2, double y2,
        double dx_y2, double x);
/*EMM-Mesh interpolation program functions*/

void EMM_cart2sphere(double x, double y, double z, double *phi_rad, double *lat_rad, double *r)
{
    double d;

    *r = sqrt(x * x + y * y + z * z);
    d = z / *r;
    *lat_rad = asin(d);
    if(x > 0)
        *phi_rad = atan(y / x);
    if(x < 0 && y >= 0)
        *phi_rad = atan(y / x) + M_PI;
    if(x < 0 && y < 0)
        *phi_rad = atan(y / x) - M_PI;
    if(fabs(x) < 1.0e-10 && y >= 0)
        *phi_rad = M_PI / 2;
    if(fabs(x) < 1.0e-10 && y < 0)
        *phi_rad = M_PI / -2;
    return;
}

void EMM_sphere2cart(double phi_rad, double lat_rad, double r, double *x, double *y, double *z)
{
    *x = r * cos(phi_rad) * cos(lat_rad);
    *y = r * sin(phi_rad) * cos(lat_rad);
    *z = r * sin(lat_rad);
    return;
}

void EMM_sphere2cart_vec(double phi_rad, double lat_rad, double vlat, double vphi, double vminusr,
        double *vx, double *vy, double *vz)
{
    double sp, st, cp, ct;

    if(fabs(fabs(lat_rad) - 0.5 * M_PI) < MIN_POLE_DIST)
    {
        vlat = 0; /* not defined at poles*/
        vphi = 0; /* not defined at poles*/
    }

    sp = sin(phi_rad);
    cp = cos(phi_rad);
    st = sin(0.5 * M_PI - lat_rad);
    ct = cos(0.5 * M_PI - lat_rad);
    *vx = -sp * vphi - ct * cp * vlat - st * cp*vminusr;
    *vy = cp * vphi - ct * sp * vlat - st * sp*vminusr;
    *vz = st * vlat - ct*vminusr;
} /* sphere2cart_vec */

void EMM_cart2sphere_vec(double phi_rad, double lat_rad, double vx, double vy, double vz,
        double *vlat, double *vphi, double *vminusr)
{
    double sp, st, cp, ct;

    sp = sin(phi_rad);
    cp = cos(phi_rad);
    st = sin(0.5 * M_PI - lat_rad);
    ct = cos(0.5 * M_PI - lat_rad);

    *vlat = -ct * cp * vx - ct * sp * vy + st*vz;
    *vphi = -sp * vx + cp*vy;
    *vminusr = -st * cp * vx - st * sp * vy - ct*vz;

    if(fabs(fabs(lat_rad) - 0.5 * M_PI) < MIN_POLE_DIST)
    {
        *vlat = 0; /* not defined at poles*/
        *vphi = 0; /* not defined at poles*/
    }
} /* cart2sphere_vec */

void EMM_geodetic2geocentric(double phi, double h, double *latrad, double *r) /* phi is geodetic latitude in radian*/
{
    double cosphi2, sinphi2, A2, B2, c;

    A2 = WGS84_A*WGS84_A;
    B2 = WGS84_B*WGS84_B;
    cosphi2 = pow(cos(phi), 2);
    sinphi2 = pow(sin(phi), 2);

    c = h * sqrt(A2 * cosphi2 + B2 * sinphi2);

    *r = sqrt(h * h + 2 * c + (A2 * A2 * cosphi2 + B2 * B2 * sinphi2) / (A2 * cosphi2 + B2 * sinphi2));

    if(fabs(phi - M_PI / 2.0) < 0.00001) *latrad = phi;
    else
    {
        *latrad = atan((c + B2) / (c + A2) * tan(phi));
    }
} /* geodetic2geocentric */

void EMM_geocentric2geodetic_vec(double theta, double delta, double bx, double by, double bz,
        double *x, double *y, double *z)
{
    double psi, sp, cp;

    psi = delta - theta; /* compatible with WMM2005 report, sign reversed for colatitudes */
    sp = sin(psi);
    cp = cos(psi);
    *x = bx * cp - bz*sp;
    *y = by;
    *z = bx * sp + bz*cp;
} /* geocentric2geodetic_vec */

static int EMM_find_alt_index(int verbose, EMM_tmesh mesh, double alt, int *ialt0, int *ialt1, double *alt0, double *alt1, double *altfac0, double *altfac1)
{

    if(alt <= (mesh.alt)[0]) /* below mesh */
    {
        *ialt0 = 0;
        *ialt1 = 1;
    } else
    {
        if(alt >= (mesh.alt)[mesh.nalt - 1]) /* above mesh */
        {
            *ialt0 = mesh.nalt - 2;
            *ialt1 = mesh.nalt - 1;
        } else /* normal case */
        {
            *ialt1 = 1;
            while(alt > (mesh.alt)[*ialt1] && *ialt1 < mesh.nalt - 1) (*ialt1)++;
            *ialt0 = *ialt1 - 1;
        }
    }

    *alt0 = (mesh.alt)[*ialt0];
    *alt1 = (mesh.alt)[*ialt1];
    *altfac0 = (*alt1 - alt) / (*alt1 - *alt0);
    *altfac1 = 1.0 - *altfac0;

    if(*altfac0 < -1.0 || *altfac0 > 1.25 || *altfac1 < -0.25 || *altfac1 > 2.0)
    {
        if(verbose)
        {
            printf("Error - Altitude %f is outside of permissible range [%4.0f to %4.0f]\n",
                    alt, (mesh.alt)[0] - 0.25 * ((mesh.alt)[1]-(mesh.alt)[0]), (mesh.alt)[mesh.nalt - 1]+((mesh.alt)[mesh.nalt - 1]-(mesh.alt)[mesh.nalt - 2]));
            fflush(stdout);
        }
        return EXIT_MESH_FILE_ALT_ERROR;
    }
    return 0;
}

static void EMM_find_lat_index(EMM_tmesh mesh, int ialt, double lat, int *ilat0, int *ilat1, double *lat0, double *lat1, double *latfac0, double *latfac1)
{

    if(lat <= -90.0 + 0.5 * (mesh.latres)[ialt]) /* South Pole */
    {
        *ilat0 = 0;
        *ilat1 = 1;
    } else
    {
        if(lat >= 90.0 - 0.5 * (mesh.latres)[ialt]) /* North Pole */
        {
            *ilat0 = (mesh.nlat)[ialt] - 2;
            *ilat1 = (mesh.nlat)[ialt] - 1;
        } else /* normal case */
        {
            *ilat0 = (int) floor((90 + lat) / (mesh.latres)[ialt] - 0.5);
            *ilat1 = *ilat0 + 1;
        }
    } /* not South Pole */

    *lat0 = (*ilat0 + 0.5) * (mesh.latres)[ialt] - 90.0;
    *lat1 = (*ilat1 + 0.5) * (mesh.latres)[ialt] - 90.0;
    *latfac0 = (*lat1 - lat) / (*lat1 - *lat0);
    *latfac1 = 1.0 - *latfac0;

    if(*latfac0 < -0.5 || *latfac0 > 1.5 || *latfac1 < -0.5 || *latfac1 > 1.5) /* allow extrapolation at poles */
    {
        printf("ERROR: latfac0=%.4f, latfac1=%.4f\n", *latfac0, *latfac1);
        exit(EXIT_FAILURE);
    }

    return;
}

static void EMM_find_lon_index(EMM_tmesh mesh, int ialt, int ilat, double lon, int *ilon0, int *ilon1, double *lon0, double *lon1, double *lonfac0, double *lonfac1)
{

    *ilon0 = (int) floor(lon / ((mesh.lonres)[ialt])[ilat]);
    *ilon1 = *ilon0 + 1;

    *lon0 = (*ilon0 * ((mesh.lonres)[ialt])[ilat]);
    *lon1 = (*ilon1 * ((mesh.lonres)[ialt])[ilat]);

    *lonfac0 = (*lon1 - lon) / (*lon1 - *lon0);
    *lonfac1 = 1.0 - *lonfac0;


    /*if (*lonfac0 < 0.0 || *lonfac0 > 1.0 || *lonfac1 < 0.0 || *lonfac1 > 1.0)
      {
        printf("ERROR: lonfac0=%.4f, lonfac1=%.4f\n",*lonfac0,*lonfac1);
        exit(EXIT_FAILURE);
      }*/

    return;
}

int EMM_mesh_convert(int verbose, char infname[], char outfname[]) /* convert a mesh file from text file format to binary format for the local hardware */
{
    int nalt, nlat, ialt, ilat, ilon, ideriv, iline = 0;
    int nlon;
    int icomp, check;
    MESH_REAL_TYPE oneval, onealt, epoch, version;
    FILE *infile, *outfile;

    if(verbose)
    {
        printf("Converting mesh file %s to binary file %s\n", infname, outfname);
        fflush(stdout);
    }

    /* open input meshfile */

    infile = fopen(infname, "r");
    if(infile == NULL)
    {
        if(verbose)
            printf("Input mesh file %s in text format could not be opened\n", infname);
        return EXIT_MESH_TEXT_FILE_NOT_FOUND;
    }

    /* open output meshfile */

    outfile = fopen(outfname, "wb");
    if(outfile == NULL)
    {
        if(verbose)
            printf("Binary output mesh file %s could not be opened\n", outfname);
        return EXIT_MESH_BIN_FILE_WRITE_ERROR;
    }

    /* read and write version and epoch */

    fscanf(infile, "%g%*[^\n]", &version);
    iline++;
    fwrite(&version, sizeof (MESH_REAL_TYPE), 1, outfile);

    fscanf(infile, "%g%*[^\n]", &epoch);
    iline++;
    fwrite(&epoch, sizeof (MESH_REAL_TYPE), 1, outfile);

    if(epoch < 1900 || epoch > 2100)
    {
        if(verbose)
        {
            printf("Error in file line %1d - Epoch out of range: %.1f \n", iline, epoch);
            fflush(stdout);
        }
        return EXIT_MESH_FILE_EPOCH_ERROR;
    }

    /* read and write number of altitude layers */

    fscanf(infile, "%d%*[^\n]", &nalt);
    iline++;
    fwrite(&nalt, sizeof (int), 1, outfile);

    if(nalt < 1 || nalt > 1000)
    {
        if(verbose)
        {
            printf("Error in file line %1d - Number of altitude layers out of range: %d \n", iline, nalt);
            fflush(stdout);
        }
        return EXIT_MESH_FILE_NALT_ERROR;
    }

    for(ialt = 0; ialt < nalt; ialt++)
    {
        fscanf(infile, "%g%*[^\n]", &onealt);
        iline++;
        fwrite(&onealt, sizeof (MESH_REAL_TYPE), 1, outfile);

        fscanf(infile, "%d%*[^\n]", &nlat);
        iline++;
        if(nlat < 1 || nlat > MESH_MAXLAT)
        {
            if(verbose)
            {
                printf("Error in file line %1d - Number of rows in latitude out of range: %d \n", iline, nlat);
                fflush(stdout);
            }
            return EXIT_MESH_FILE_NLAT_ERROR;
        }
        fwrite(&nlat, sizeof (int), 1, outfile);

        for(ilat = 0; ilat < nlat; ilat++)
        {
            fscanf(infile, "%d%*[^\n]", &nlon);
            iline++;
            if(nlon < 1 || nlon > 2 * MESH_MAXLAT)
            {
                if(verbose)
                {
                    printf("Error in file line %1d - Number of cells at altitude %1d in row %1d is out of range: nlon[%1d][%1d]=%1d \n",
                            iline, ialt, ilat, ialt, ilat, nlon);
                    fflush(stdout);
                }
                return EXIT_MESH_FILE_NLON_ERROR;
            }

            fwrite(&nlon, sizeof (int), 1, outfile);

            /* copying data */

            for(ilon = 0; ilon < nlon; ilon++)
            {
                for(icomp = 0; icomp < 3; icomp++)
                {
                    for(ideriv = 0; ideriv < 4; ideriv++)
                    {
                        fscanf(infile, "%g", &oneval);
                        fwrite(&oneval, sizeof (MESH_REAL_TYPE), 1, outfile);
                    }
                    fscanf(infile, "%*[^\n]");
                    iline++;
                } /* for icomp */
            } /* for ilon */
        } /* for ilat */
    } /* for ialt */

    check = 4711;
    fwrite(&check, sizeof (int), 1, outfile);

    fclose(infile);
    fclose(outfile);

    if(verbose)
    {
        printf("Completed converting mesh file to binary\n");
        fflush(stdout);
    }
    return 0;
} /* mesh_convert */

int EMM_mesh_read(int verbose, char meshfname[], EMM_tmesh *mesh) /* read binary mesh file, allocate memory and fill mesh */
{
    int ialt, ilat, ilon;
    int icomp, ideriv, check;
    MESH_REAL_TYPE oneval;
    FILE *meshfile;


    /* open meshfile */

    if(verbose)
    {
        printf("Reading binary mesh file %s...\n", meshfname);
        fflush(stdout);
    }
    meshfile = fopen(meshfname, "rb");
    if(meshfile == NULL)
    {
        if(verbose)
            printf("Binary mesh file %s could not be opened\n", meshfname);
        return EXIT_MESH_BIN_FILE_NOT_FOUND;
    }

    /* read version and epoch */

    fread(&oneval, sizeof (MESH_REAL_TYPE), 1, meshfile);
    ((*mesh).version) = oneval;
    fread(&oneval, sizeof (MESH_REAL_TYPE), 1, meshfile);
    ((*mesh).epoch) = oneval;
    if((*mesh).epoch < 1900 || (*mesh).epoch > 2100)
    {
        if(verbose)
        {
            printf("Error - Epoch out of range: %.1f \n", (*mesh).epoch);
            fflush(stdout);
        }
        return EXIT_MESH_FILE_EPOCH_ERROR;
    }

    /* read number of altitude layers */

    fread(&((*mesh).nalt), sizeof (int), 1, meshfile);
    if((*mesh).nalt < 1 || (*mesh).nalt > 1000)
    {
        if(verbose)
        {
            printf("Error - Number of altitude layers out of range: %d \n", (*mesh).nalt);
            fflush(stdout);
        }
        return EXIT_MESH_FILE_NALT_ERROR;
    }

    /* allocate space for altitudes */

    (*mesh).alt = malloc((*mesh).nalt * sizeof (double));
    if((*mesh).alt == NULL)
    {
        if(verbose)
        {
            printf("Error - out of memory\n");
            fflush(stdout);
        }
        return EXIT_MESH_MEM_ALLOC_ERROR;
    }

    /* allocate space for number of latitudes and latitude resolution in each layer */

    (*mesh).nlat = malloc((*mesh).nalt * sizeof (int));
    if((*mesh).nlat == NULL)
    {
        if(verbose)
        {
            printf("Error - out of memory\n");
            fflush(stdout);
        }
        return EXIT_MESH_MEM_ALLOC_ERROR;
    }

    (*mesh).latres = malloc((*mesh).nalt * sizeof (double));
    if((*mesh).latres == NULL)
    {
        if(verbose)
        {
            printf("Error - out of memory\n");
            fflush(stdout);
        }
        return EXIT_MESH_MEM_ALLOC_ERROR;
    }

    /* allocate space for pointer to number of longitudes and longitude resolution for each latitude in each layer */

    (*mesh).nlon = malloc((*mesh).nalt * sizeof (int*));
    if((*mesh).nlon == NULL)
    {
        if(verbose)
        {
            printf("Error - out of memory\n");
            fflush(stdout);
        }
        return EXIT_MESH_MEM_ALLOC_ERROR;
    }

    (*mesh).lonres = malloc((*mesh).nalt * sizeof (double**));
    if((*mesh).lonres == NULL)
    {
        if(verbose)
        {
            printf("Error - out of memory\n");
            fflush(stdout);
        }
        return EXIT_MESH_MEM_ALLOC_ERROR;
    }

    /* allocate data pointers for each altitude layer */

    (*mesh).comp = malloc((*mesh).nalt * sizeof (double***));
    if((*mesh).comp == NULL)
    {
        if(verbose)
        {
            printf("Error - out of memory\n");
            fflush(stdout);
        }
        return EXIT_MESH_MEM_ALLOC_ERROR;
    }

    for(ialt = 0; ialt < (*mesh).nalt; ialt++)
    {

        /* get layer altitude */

        fread(&oneval, sizeof (MESH_REAL_TYPE), 1, meshfile);
        ((*mesh).alt)[ialt] = oneval;

        /* get number of latitudes for this layer */
        /* different altitudes can have different resolution */

        fread(&((*mesh).nlat)[ialt], sizeof (int), 1, meshfile);
        if(((*mesh).nlat)[ialt] < 1 || ((*mesh).nlat)[ialt] > MESH_MAXLAT)
        {
            if(verbose)
            {
                printf("Error - Number of rows in latitude out of range: %d \n", ((*mesh).nlat)[ialt]);
                fflush(stdout);
            }
            return EXIT_MESH_FILE_NLAT_ERROR;
        }

        ((*mesh).latres)[ialt] = 180.0 / ((*mesh).nlat)[ialt];

        /* allocate space for number of longitudes and resolution for each latitude in each layer */

        ((*mesh).nlon)[ialt] = malloc(((*mesh).nlat)[ialt] * sizeof (int*));
        if(((*mesh).nlon)[ialt] == NULL)
        {
            if(verbose)
            {
                printf("Error - out of memory\n");
                fflush(stdout);
            }
            return EXIT_MESH_MEM_ALLOC_ERROR;
        }

        ((*mesh).lonres)[ialt] = malloc(((*mesh).nlat)[ialt] * sizeof (double*));
        if(((*mesh).lonres)[ialt] == NULL)
        {
            if(verbose)
            {
                printf("Error - out of memory\n");
                fflush(stdout);
            }
            return EXIT_MESH_MEM_ALLOC_ERROR;
        }


        /* allocate data pointer for each latitude at this altitude */

        ((*mesh).comp)[ialt] = malloc(((*mesh).nlat)[ialt] * sizeof (double**));
        if(((*mesh).comp)[ialt] == NULL)
        {
            if(verbose)
            {
                printf("Error - out of memory\n");
                fflush(stdout);
            }
            return EXIT_MESH_MEM_ALLOC_ERROR;
        }

        for(ilat = 0; ilat < ((*mesh).nlat)[ialt]; ilat++)
        {
            fread(&(((*mesh).nlon)[ialt])[ilat], sizeof (int), 1, meshfile);
            if((((*mesh).nlon)[ialt])[ilat] < 1 || (((*mesh).nlon)[ialt])[ilat] > 2 * ((*mesh).nlat)[ialt])
            {
                if(verbose)
                {
                    printf("Error - Number of cells at altitude %1d in row %1d is out of range: nlon[%1d][%1d]=%1d \n", ialt, ilat, ialt, ilat, (((*mesh).nlon)[ialt])[ilat]);
                    fflush(stdout);
                }
                return EXIT_MESH_FILE_NLON_ERROR;
            }
            (((*mesh).lonres)[ialt])[ilat] = 360.0 / (((*mesh).nlon)[ialt])[ilat];

            /* allocate data pointer for each latitude at this altitude */

            (((*mesh).comp)[ialt])[ilat] = malloc(((((*mesh).nlon)[ialt])[ilat] + 1) * sizeof (double*)); /* add one to duplicate the first longitude at the end */

            if((((*mesh).comp)[ialt])[ilat] == NULL)
            {
                if(verbose)
                {
                    printf("Error - out of memory\n");
                    fflush(stdout);
                }
                return EXIT_MESH_MEM_ALLOC_ERROR;
            } /* if NULL */

            for(ilon = 0; ilon < (((*mesh).nlon)[ialt])[ilat]; ilon++)
            {
                ((((*mesh).comp)[ialt])[ilat])[ilon] = malloc(12 * sizeof (double));
                if(((((*mesh).comp)[ialt])[ilat])[ilon] == NULL)
                {
                    if(verbose)
                    {
                        printf("Error - out of memory\n");
                        fflush(stdout);
                    }
                    return EXIT_MESH_MEM_ALLOC_ERROR;
                } /* if NULL */


                for(icomp = 0; icomp < 3; icomp++)
                    for(ideriv = 0; ideriv < 4; ideriv++)
                    {
                        fread(&oneval, sizeof (MESH_REAL_TYPE), 1, meshfile);
                        (((((*mesh).comp)[ialt])[ilat])[ilon])[icomp * 4 + ideriv] = oneval;
                    }
            } /* for ilon */

            /* duplicate the first longitude at the end */

            ilon = (((*mesh).nlon)[ialt])[ilat]; /* now ilon = nlon for that row */
            ((((*mesh).comp)[ialt])[ilat])[ilon] = malloc(12 * sizeof (double)); /* add one to duplicate the first longitude at the end */
            if(((((*mesh).comp)[ialt])[ilat])[ilon] == NULL)
            {
                if(verbose)
                {
                    printf("Error - out of memory\n");
                    fflush(stdout);
                }
                return EXIT_MESH_MEM_ALLOC_ERROR;
            } /* if NULL */

            for(icomp = 0; icomp < 3; icomp++)
                for(ideriv = 0; ideriv < 4; ideriv++)
                    (((((*mesh).comp)[ialt])[ilat])[ilon])[icomp * 4 + ideriv] = (((((*mesh).comp)[ialt])[ilat])[0])[icomp * 4 + ideriv];

        } /* for ilat */
    } /* for ialt */

    fread(&check, sizeof (int), 1, meshfile);
    if(check != 4711)
    {
        if(verbose)
        {
            printf("Mesh file format error. Last integer in binary file should be 4711, but is = %1d\n", check);
            fflush(stdout);
        }
        return EXIT_MESH_FILE_FORMAT_ERROR;
    }


    if(verbose) printf("Mesh has %1d altitude layers\n", (*mesh).nalt);

    fclose(meshfile);

    if(verbose)
    {
        printf("Completed reading mesh file\n");
        fflush(stdout);
    }

    return 0;
} /* mesh_read */

static double EMM_cell(EMM_tmesh mesh, int ialt, int ilat, int ilon, int icomp)
{
    return ((((mesh.comp)[ialt])[ilat])[ilon])[icomp * 4];
}

static double EMM_dlon_cell(EMM_tmesh mesh, int ialt, int ilat, int ilon, int icomp)
{
    return ((((mesh.comp)[ialt])[ilat])[ilon])[icomp * 4 + 1];
}

static double EMM_dlat_cell(EMM_tmesh mesh, int ialt, int ilat, int ilon, int icomp)
{
    return ((((mesh.comp)[ialt])[ilat])[ilon])[icomp * 4 + 2];
}

static double EMM_dalt_cell(EMM_tmesh mesh, int ialt, int ilat, int ilon, int icomp)
{
    return ((((mesh.comp)[ialt])[ilat])[ilon])[icomp * 4 + 3];
}

static double EMM_interpolate_quadratic(double x1, double y1, double dx_y1, double x2, double y2, double dx_y2, double x)
/* interpolating function: y = ax^2 + bx + c */
{
    double a, b, c, dx;

    dx = x2 - x1;
    a = (dx_y2 - dx_y1) / (2.0 * dx);
    b = (y2 - y1 - a * (x2 * x2 - x1 * x1)) / dx;
    c = y1 - a * x1 * x1 - b*x1;
    return (a * x * x + b * x + c);
} /* interpolate_quadratic */

static double EMM_interpolate_cubic(double x1, double y1, double dx_y1, double x2, double y2, double dx_y2, double x)
/* interpolating function: y = ax^3 + bx^2 + cx + d */
{
    double a, b, c, d, dx, dy, dydx, x1s, x1c, x2s, x2c;

    dx = x2 - x1;
    dy = y2 - y1;
    dydx = dy / dx;
    x1s = x1*x1;
    x1c = x1s*x1;
    x2s = x2*x2;
    x2c = x2s*x2;

    a = (dx_y1 + dx_y2 - 2.0 * dydx) / (dx * dx);

    b = (dx_y1 - dydx - (2.0 * x1s - x1 * x2 - x2s) * a) / (-dx);

    c = (dy - (x2c - x1c) * a - (x2s - x1s) * b) / dx;

    d = y1 - x1c * a - x1s * b - x1*c;

    return (a * x * x * x + b * x * x + c * x + d);
} /* interpolate_cubic */

static double EMM_interpolate(double x1, double y1, double dx_y1, double x2, double y2, double dx_y2, double x)
{
    if(INTERPOLATE_CUBIC)
        return EMM_interpolate_cubic(x1, y1, dx_y1, x2, y2, dx_y2, x);
    else
        return EMM_interpolate_quadratic(x1, y1, dx_y1, x2, y2, dx_y2, x);
} /* interpolate */

int EMM_mesh_interpolate(int verbose, EMM_tmesh mesh, double lon, double lat, double alt, double geoc_lat, double *geoc_Bx, double *geoc_By, double *geoc_Bz)
{
    int icomp;
    int ialt0, ialt1;
    int ilat00, ilat01, ilat10, ilat11;
    int ilon000, ilon001, ilon010, ilon011, ilon100, ilon101, ilon110, ilon111;

    double alt0, alt1;
    double lat00, lat01, lat10, lat11;
    double lon000, lon001, lon010, lon011, lon100, lon101, lon110, lon111;

    double altfac0, altfac1;
    double latfac00, latfac01, latfac10, latfac11;
    double lonfac000, lonfac001, lonfac010, lonfac011, lonfac100, lonfac101, lonfac110, lonfac111;

    double dlat_cell0, dlat_cell1;
    double dalt_cell00, dalt_cell01, dalt_cell10, dalt_cell11, dalt_cell0, dalt_cell1;

    double f0, f1, g0, g1, h[3];
    double lon_rad, geoc_lat_rad;

    while(lon < 0) lon += 360;
    while(lon >= 360) lon -= 360;

    EMM_find_alt_index(verbose, mesh, alt, &ialt0, &ialt1, &alt0, &alt1, &altfac0, &altfac1);

    /* lower altitude */

    EMM_find_lat_index(mesh, ialt0, lat, &ilat00, &ilat01, &lat00, &lat01, &latfac00, &latfac01);

    /* lower altitude, lower lat */

    EMM_find_lon_index(mesh, ialt0, ilat00, lon, &ilon000, &ilon001, &lon000, &lon001, &lonfac000, &lonfac001);

    /* lower altitude, upper lat */

    EMM_find_lon_index(mesh, ialt0, ilat01, lon, &ilon010, &ilon011, &lon010, &lon011, &lonfac010, &lonfac011);

    /* upper altitude */

    EMM_find_lat_index(mesh, ialt1, lat, &ilat10, &ilat11, &lat10, &lat11, &latfac10, &latfac11);

    /* upper altitude, lower lat */

    EMM_find_lon_index(mesh, ialt1, ilat10, lon, &ilon100, &ilon101, &lon100, &lon101, &lonfac100, &lonfac101);

    /* upper altitude, upper lat */

    EMM_find_lon_index(mesh, ialt1, ilat11, lon, &ilon110, &ilon111, &lon110, &lon111, &lonfac110, &lonfac111);

    /* interpolate */

    for(icomp = 0; icomp < 3; icomp++)
    {

#if USE_DERIVATIVES

        /* lower altitude */

        f0 = EMM_interpolate(lon000, EMM_cell(mesh, ialt0, ilat00, ilon000, icomp), EMM_dlon_cell(mesh, ialt0, ilat00, ilon000, icomp),
                lon001, EMM_cell(mesh, ialt0, ilat00, ilon001, icomp), EMM_dlon_cell(mesh, ialt0, ilat00, ilon001, icomp), lon);

        f1 = EMM_interpolate(lon010, EMM_cell(mesh, ialt0, ilat01, ilon010, icomp), EMM_dlon_cell(mesh, ialt0, ilat01, ilon010, icomp),
                lon011, EMM_cell(mesh, ialt0, ilat01, ilon011, icomp), EMM_dlon_cell(mesh, ialt0, ilat01, ilon011, icomp), lon);

        dlat_cell0 = lonfac000 * EMM_dlat_cell(mesh, ialt0, ilat00, ilon000, icomp) + lonfac001 * EMM_dlat_cell(mesh, ialt0, ilat00, ilon001, icomp);
        dlat_cell1 = lonfac010 * EMM_dlat_cell(mesh, ialt0, ilat01, ilon010, icomp) + lonfac011 * EMM_dlat_cell(mesh, ialt0, ilat01, ilon011, icomp);

        g0 = EMM_interpolate(lat00, f0, dlat_cell0, lat01, f1, dlat_cell1, lat);

        dalt_cell00 = lonfac000 * EMM_dalt_cell(mesh, ialt0, ilat00, ilon000, icomp) + lonfac001 * EMM_dalt_cell(mesh, ialt0, ilat00, ilon001, icomp);
        dalt_cell01 = lonfac010 * EMM_dalt_cell(mesh, ialt0, ilat01, ilon010, icomp) + lonfac011 * EMM_dalt_cell(mesh, ialt0, ilat01, ilon011, icomp);
        dalt_cell0 = latfac00 * dalt_cell00 + latfac01*dalt_cell01;

        /* upper altitude */

        f0 = EMM_interpolate(lon100, EMM_cell(mesh, ialt1, ilat10, ilon100, icomp), EMM_dlon_cell(mesh, ialt1, ilat10, ilon100, icomp),
                lon101, EMM_cell(mesh, ialt1, ilat10, ilon101, icomp), EMM_dlon_cell(mesh, ialt1, ilat10, ilon101, icomp), lon);

        f1 = EMM_interpolate(lon110, EMM_cell(mesh, ialt1, ilat11, ilon110, icomp), EMM_dlon_cell(mesh, ialt1, ilat11, ilon110, icomp),
                lon111, EMM_cell(mesh, ialt1, ilat11, ilon111, icomp), EMM_dlon_cell(mesh, ialt1, ilat11, ilon111, icomp), lon);

        dlat_cell0 = lonfac100 * EMM_dlat_cell(mesh, ialt1, ilat10, ilon100, icomp) + lonfac101 * EMM_dlat_cell(mesh, ialt1, ilat10, ilon101, icomp);
        dlat_cell1 = lonfac110 * EMM_dlat_cell(mesh, ialt1, ilat11, ilon110, icomp) + lonfac111 * EMM_dlat_cell(mesh, ialt1, ilat11, ilon111, icomp);

        g1 = EMM_interpolate(lat10, f0, dlat_cell0, lat11, f1, dlat_cell1, lat);

        dalt_cell10 = lonfac100 * EMM_dalt_cell(mesh, ialt1, ilat10, ilon100, icomp) + lonfac101 * EMM_dalt_cell(mesh, ialt1, ilat10, ilon101, icomp);
        dalt_cell11 = lonfac110 * EMM_dalt_cell(mesh, ialt1, ilat11, ilon110, icomp) + lonfac111 * EMM_dalt_cell(mesh, ialt1, ilat11, ilon111, icomp);
        dalt_cell1 = latfac10 * dalt_cell10 + latfac11*dalt_cell11;

        /* final value */

        h[icomp] = EMM_interpolate(alt0, g0, dalt_cell0, alt1, g1, dalt_cell1, alt);

#else  /* linear interpolation of values without taking curvature into account */

        f0 = lonfac000 * EMM_cell(mesh, ialt0, ilat00, ilon000, icomp) + lonfac001 * EMM_cell(mesh, ialt0, ilat00, ilon001, icomp);
        f1 = lonfac010 * EMM_cell(mesh, ialt0, ilat01, ilon010, icomp) + lonfac011 * EMM_cell(mesh, ialt0, ilat01, ilon011, icomp);
        g0 = latfac00 * f0 + latfac01*f1;

        f0 = lonfac100 * EMM_cell(mesh, ialt1, ilat10, ilon100, icomp) + lonfac101 * EMM_cell(mesh, ialt1, ilat10, ilon101, icomp);
        f1 = lonfac110 * EMM_cell(mesh, ialt1, ilat11, ilon110, icomp) + lonfac111 * EMM_cell(mesh, ialt1, ilat11, ilon111, icomp);
        g1 = latfac10 * f0 + latfac11*f1;

        h[icomp] = altfac0 * g0 + altfac1*g1;

#endif   /* if linear interpolation */

    } /* for icomp */


    /* h now contains the interpolated vector components in a global cartesian reference frame */
    /* these components now have to be transformed to the local geocentric frame. Later, the user */
    /* still has to transform them further into the local geographic (WGS84) frame  */


    /* Transform the cartesian vector components into the local geocentric frame */

    lon_rad = lon * M_PI / 180.0;
    geoc_lat_rad = geoc_lat * M_PI / 180.0;

    EMM_cart2sphere_vec(lon_rad, geoc_lat_rad, h[0], h[1], h[2], geoc_Bx, geoc_By, geoc_Bz); /* x=north, y=east, z=down */

    return 0;
}

int EMM_PointCalcFromMesh(MAGtype_CoordGeodetic CordGeo, MAGtype_CoordSpherical CordSph, MAGtype_Date UserDate, MAGtype_MagneticResults *MagResults, EMM_tmesh mesh, EMM_tmesh mesh_SV)
{
    int verbose = 1;
    double lon, lon_rad, geod_lat, geod_lat_rad, alt, date, geoc_lat, r, x, y, z, geod_colat_rad, geoc_colat_rad, xs, ys, zs, geod_x, geod_y, geod_z;

    lon = CordGeo.lambda;
    lon_rad = lon * M_PI / 180.0;

    geod_lat = CordGeo.phi;
    geod_lat_rad = geod_lat * M_PI / 180.0;

    alt = CordGeo.HeightAboveEllipsoid;

    date = UserDate.DecimalYear;

    geoc_lat = CordSph.phig;

    r = CordSph.r;



    EMM_mesh_interpolate(verbose, mesh, lon, geod_lat, alt, geoc_lat, &x, &y, &z);
    EMM_mesh_interpolate(verbose, mesh_SV, lon, geod_lat, alt, geoc_lat, &xs, &ys, &zs);

    x += (date - mesh.epoch) * xs;
    y += (date - mesh.epoch) * ys;
    z += (date - mesh.epoch) * zs;

    geod_colat_rad = (90.0 - geod_lat) * M_PI / 180.0;
    geoc_colat_rad = (90.0 - geoc_lat) * M_PI / 180.0;
    EMM_geocentric2geodetic_vec(geoc_colat_rad, geod_colat_rad, x, y, z, &geod_x, &geod_y, &geod_z);

    MagResults->Bx = geod_x;
    MagResults->By = geod_y;
    MagResults->Bz = geod_z;

    return 1;
} /*EMM_PointCalcFromMesh*/

int EMM_Grid(MAGtype_CoordGeodetic minimum, MAGtype_CoordGeodetic maximum, double
        cord_step_size, double altitude_step_size, double time_step, MAGtype_Geoid
        *Geoid, MAGtype_Ellipsoid Ellip, MAGtype_Date StartDate, MAGtype_Date EndDate, int ElementOption, int PrintOption, char *OutputFile, EMM_tmesh mesh, EMM_tmesh mesh_SV)
{
    double a, b, c, d, PrintElement;

    MAGtype_CoordSpherical CoordSpherical;
    MAGtype_MagneticResults MagneticResultsGeo, MagneticVariation;
    MAGtype_GeoMagneticElements GeoMagneticElements;


    FILE *fileout;

    if(PrintOption == 1)
    {
        fileout = fopen(OutputFile, "w");
        if(!fileout)
        {
            printf("Error opening %s to write", OutputFile);
            return FALSE;
        }
    }


    if(fabs(cord_step_size) < 1.0e-10) cord_step_size = 99999.0; //checks to make sure that the step_size is not too small
    if(fabs(altitude_step_size) < 1.0e-10) altitude_step_size = 99999.0;
    if(fabs(time_step) < 1.0e-10) time_step = 99999.0;


    a = minimum.HeightAboveGeoid; //sets the loop intialization values
    b = minimum.phi;
    c = minimum.lambda;
    d = StartDate.DecimalYear;



    for(minimum.HeightAboveGeoid = a; minimum.HeightAboveGeoid <= maximum.HeightAboveGeoid; minimum.HeightAboveGeoid += altitude_step_size) /* Altitude loop*/
    {


        for(minimum.phi = b; minimum.phi <= maximum.phi; minimum.phi += cord_step_size) /*Latitude loop*/
        {



            for(minimum.lambda = c; minimum.lambda <= maximum.lambda; minimum.lambda += cord_step_size) /*Longitude loop*/
            {
                if(Geoid->UseGeoid == 1)
                    MAG_ConvertGeoidToEllipsoidHeight(&minimum, Geoid); //This converts the height above mean sea level to height above the WGS-84 ellipsoid
                else
                    minimum.HeightAboveEllipsoid = minimum.HeightAboveGeoid;
                MAG_GeodeticToSpherical(Ellip, minimum, &CoordSpherical);
                for(StartDate.DecimalYear = d; StartDate.DecimalYear <= EndDate.DecimalYear; StartDate.DecimalYear += time_step) /*Year loop*/
                {

                    EMM_PointCalcFromMesh(minimum, CoordSpherical, StartDate, &MagneticResultsGeo, mesh, mesh_SV);
                    MAG_CalculateGeoMagneticElements(&MagneticResultsGeo, &GeoMagneticElements); /* Calculate the Geomagnetic elements, Equation 18 , WMM Technical report */
                    /*For calculating secular variation*/
                    StartDate.DecimalYear += 1;
                    EMM_PointCalcFromMesh(minimum, CoordSpherical, StartDate, &MagneticVariation, mesh, mesh_SV);
                    MagneticVariation.Bx += -MagneticResultsGeo.Bx;
                    MagneticVariation.By += -MagneticResultsGeo.By;
                    MagneticVariation.Bz += -MagneticResultsGeo.Bz;
                    MAG_CalculateSecularVariationElements(MagneticVariation, &GeoMagneticElements); /*Calculate the secular variation of each of the Geomagnetic elements, Equation 19, WMM Technical report*/
                    /*Return to correct year*/
                    StartDate.DecimalYear -= 1;
                    switch(ElementOption) {
                        case 1:
                            PrintElement = GeoMagneticElements.Decl; /*1. Angle between the magnetic field vector and true north, positive east*/
                            break;
                        case 2:
                            PrintElement = GeoMagneticElements.Incl; /*2. Angle between the magnetic field vector and the horizontal plane, positive downward*/
                            break;
                        case 3:
                            PrintElement = GeoMagneticElements.F; /*3. Magnetic Field Strength*/
                            break;
                        case 4:
                            PrintElement = GeoMagneticElements.H; /*4. Horizontal Magnetic Field Strength*/
                            break;
                        case 5:
                            PrintElement = GeoMagneticElements.X; /*5. Northern component of the magnetic field vector*/
                            break;
                        case 6:
                            PrintElement = GeoMagneticElements.Y; /*6. Eastern component of the magnetic field vector*/
                            break;
                        case 7:
                            PrintElement = GeoMagneticElements.Z; /*7. Downward component of the magnetic field vector*/
                            break;
                        case 8:
                            PrintElement = GeoMagneticElements.GV; /*8. The Grid Variation*/
                            break;
                        case 9:
                            PrintElement = GeoMagneticElements.Decldot; /*9. Yearly Rate of change in declination*/
                            break;
                        case 10:
                            PrintElement = GeoMagneticElements.Incldot; /*10. Yearly Rate of change in inclination*/
                            break;
                        case 11:
                            PrintElement = GeoMagneticElements.Fdot; /*11. Yearly rate of change in Magnetic field strength*/
                            break;
                        case 12:
                            PrintElement = GeoMagneticElements.Hdot; /*12. Yearly rate of change in horizontal field strength*/
                            break;
                        case 13:
                            PrintElement = GeoMagneticElements.Xdot; /*13. Yearly rate of change in the northern component*/
                            break;
                        case 14:
                            PrintElement = GeoMagneticElements.Ydot; /*14. Yearly rate of change in the eastern component*/
                            break;
                        case 15:
                            PrintElement = GeoMagneticElements.Zdot; /*15. Yearly rate of change in the downward component*/
                            break;
                        case 16:
                            PrintElement = GeoMagneticElements.GVdot;
                            /*16. Yearly rate of chnage in grid variation*/;
                            break;
                        default:
                            PrintElement = GeoMagneticElements.Decl; /* 1. Angle between the magnetic field vector and true north, positive east*/
                    }

                    if(Geoid->UseGeoid == 1)
                    {
                        if(PrintOption == 1) fprintf(fileout, "%5.2lf %6.2lf %8.4lf %7.2lf %10.2lf\n", minimum.phi, minimum.lambda, minimum.HeightAboveGeoid, StartDate.DecimalYear, PrintElement);
                        else printf("%5.2lf %6.2lf %8.4lf %7.2lf %10.2lf\n", minimum.phi, minimum.lambda, minimum.HeightAboveGeoid, StartDate.DecimalYear, PrintElement);
                    } else
                    {
                        if(PrintOption == 1) fprintf(fileout, "%5.2lf %6.2lf %8.4lf %7.2lf %10.2lf\n", minimum.phi, minimum.lambda, minimum.HeightAboveEllipsoid, StartDate.DecimalYear, PrintElement);
                        else printf("%5.2lf %6.2lf %8.4lf %7.2lf %10.2lf\n", minimum.phi, minimum.lambda, minimum.HeightAboveEllipsoid, StartDate.DecimalYear, PrintElement);
                    }





                } /* year loop */

            } /*Longitude Loop */

        } /* Latitude Loop */

    } /* Altitude Loop */
    if(PrintOption == 1) fclose(fileout);



    return TRUE;
} /*EMM_Grid*/
