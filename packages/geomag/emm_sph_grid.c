
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include "GeomagnetismHeader.h"
#include "EGM9615.h"

/*
 *
 *  Revision Number: $Revision: 1324 $
 *  Last changed by: $Author: awoods $
 *  Last changed on: $Date: 2015-05-13 16:04:23 -0600 (Wed, 13 May 2015) $
 *  Manoj.C.Nair@NOAA.Gov

 */



int main()
{
    MAGtype_MagneticModel * MagneticModels[17], *TimedMagneticModel;
    MAGtype_Ellipsoid Ellip;
    MAGtype_CoordGeodetic minimum, maximum;
    MAGtype_CoordSpherical CoordSpherical;
    MAGtype_Geoid Geoid;
    MAGtype_Date startdate, enddate;
    MAGtype_LegendreFunction *LegendreFunction;
    MAGtype_SphericalHarmonicVariables *SphVariables;
    MAGtype_MagneticResults MagneticResultsSph, MagneticResultsSphVar, MagneticResultsGeo, MagneticResultsGeoVar;
    MAGtype_GeoMagneticElements GeoMagneticElements;
    int ElementOption, PrintOption, i, Epoch, nMaxEMM, NumTerms, LoadedEpoch = -1, nMax, index;
    int PrintFullGradient = 0;
    const int epochs = 16;
    double cord_step_size, altitude_step_size, time_step_size, a, b, c, d, PrintElement = 1;
    char filename[] = "EMM2015.COF";
    char filenameSV[] = "EMM2015SV.COF";
    char OutputFilename[32];
    char VersionDate_Large[] = "$Date: 2015-05-13 16:04:23 -0600 (Wed, 13 May 2015) $";
    char VersionDate[12];
    char ans[20];
    FILE *fileout;
    MAGtype_Gradient Gradient;

    strncpy(VersionDate, VersionDate_Large + 39, 11);
    VersionDate[11] = '\0';
    for(Epoch = 0; Epoch < epochs; Epoch++)
    {
        sprintf(filename, "EMM%d.COF", Epoch + 2000);
        sprintf(filenameSV, "EMM%dSV.COF", Epoch + 2000);
        if(Epoch == epochs-1)
            Epoch++;
        if(!MAG_robustReadMagneticModel_Large(filename, filenameSV, &MagneticModels[Epoch]))  {
            printf("\n EMM%d.COF or EMM%dSV.COF not found.  Press enter to exit... \n ", Epoch+2000, Epoch+2000);
            fgets(ans, 20, stdin);
            return 1;
        }
    }
    /*Create Main Field Model for EMM2015*/
    nMaxEMM = MagneticModels[0]->nMax;
    NumTerms = ((nMaxEMM + 1) * (nMaxEMM + 2) / 2);
    MagneticModels[epochs - 1] = MAG_AllocateModelMemory(NumTerms);
    MagneticModels[epochs - 1]->nMax = nMaxEMM;
    MagneticModels[epochs - 1]->nMaxSecVar = nMaxEMM;
    MagneticModels[epochs - 1]->epoch = MagneticModels[0]->epoch + epochs - 1;
    MAG_AssignMagneticModelCoeffs(MagneticModels[epochs - 1], MagneticModels[epochs], MagneticModels[epochs - 1]->nMax, MagneticModels[epochs - 1]->nMaxSecVar);

    nMax = MagneticModels[epochs]->nMax;
    NumTerms = (nMax + 1) * (nMax + 2) / 2;
    TimedMagneticModel = MAG_AllocateModelMemory(NumTerms);
    LegendreFunction = MAG_AllocateLegendreFunctionMemory(NumTerms); /* For storing the ALF functions */
    SphVariables = MAG_AllocateSphVarMemory(nMax);
    MAG_SetDefaults(&Ellip, &Geoid);
    //MAG_InitializeGeoid(&Geoid); Deprecated
    /* Set EGM96 Geoid parameters */
    Geoid.GeoidHeightBuffer = GeoidHeightBuffer;
    Geoid.Geoid_Initialized = 1;
    /* Set EGM96 Geoid parameters END */

    printf("\n\n Welcome to the Enhanced Magnetic Model (EMM) C-Program\n");
    printf("of the US National Geophysical Data Center\n");
    printf("\t\t--- Grid Calculation Program ----\n\t     --- Model Release Year: %d ---\n\t     --- Software Release Date: %s ---\n", (int)MagneticModels[epochs-1]->epoch, VersionDate);
    printf("This program may be used to generate a grid/volume of magnetic field values\n");
    printf("over latitude, longitude, altitude and time axes. To skip an axis,\n");
    printf("keep the start and end values the same and enter zero for the step size.\n");

    printf("\n\n                     Enter grid parameters \n\n");


    /* Get the Lat/Long, Altitude, Time limits from a user interface and print the grid to screen */

    MAG_GetUserGrid(&minimum, &maximum, &cord_step_size, &altitude_step_size, &time_step_size, &startdate, &enddate, &ElementOption, &PrintOption, OutputFilename, &Geoid);
    if(PrintOption == 1)
    {
        fileout = fopen(OutputFilename, "w");
        if(!fileout)
        {
            printf("Error opening %s to write", OutputFilename);
            return FALSE;
        }
    }

    if(fabs(cord_step_size) < 1.0e-10) cord_step_size = 99999.0; //checks to make sure that the step_size is not too small
    if(fabs(altitude_step_size) < 1.0e-10) altitude_step_size = 99999.0;
    if(fabs(time_step_size) < 1.0e-10) time_step_size = 99999.0;

    a = minimum.HeightAboveGeoid; //sets the loop initialization values
    b = minimum.phi;
    c = minimum.lambda;
    d = startdate.DecimalYear;
   
    for(minimum.HeightAboveGeoid = a; minimum.HeightAboveGeoid <= maximum.HeightAboveGeoid; minimum.HeightAboveGeoid += altitude_step_size) /* Altitude loop*/
    {
        for(minimum.phi = b; minimum.phi <= maximum.phi; minimum.phi += cord_step_size) /*Latitude loop*/
        {
            for(minimum.lambda = c; minimum.lambda <= maximum.lambda; minimum.lambda += cord_step_size) /*Longitude loop*/
            {
                if(Geoid.UseGeoid == 1)
                    MAG_ConvertGeoidToEllipsoidHeight(&minimum, &Geoid); //This converts the height above mean sea level to height above the WGS-84 ellipsoid
                else
                    minimum.HeightAboveEllipsoid = minimum.HeightAboveGeoid;
                MAG_GeodeticToSpherical(Ellip, minimum, &CoordSpherical);
                MAG_ComputeSphericalHarmonicVariables(Ellip, CoordSpherical, nMax, SphVariables); /* Compute Spherical Harmonic variables  */
                MAG_AssociatedLegendreFunction(CoordSpherical, nMax, LegendreFunction); /* Compute ALF  Equations 5-6, WMM Technical report*/

                for(startdate.DecimalYear = d; startdate.DecimalYear <= enddate.DecimalYear; startdate.DecimalYear += time_step_size) /*Year loop*/
                {
                    Epoch = ((int) startdate.DecimalYear - MagneticModels[0]->epoch);
                    if(Epoch < 0) Epoch = 0;
                    if(Epoch > epochs - 1) Epoch = epochs - 1;
                    if(LoadedEpoch != Epoch)
                    {
                        MagneticModels[epochs]->epoch = MagneticModels[Epoch]->epoch;
                        MAG_AssignMagneticModelCoeffs(MagneticModels[epochs], MagneticModels[Epoch], MagneticModels[Epoch]->nMax, MagneticModels[Epoch]->nMaxSecVar);
                        if(Epoch < epochs - 1)
                        {
                            for(i = 0; i < 16; i++)
                            {
                                index = 16 * 17 / 2 + i;
                                MagneticModels[epochs]->Secular_Var_Coeff_G[index] = 0;
                                MagneticModels[epochs]->Secular_Var_Coeff_H[index] = 0;
                            }
                        }
                        LoadedEpoch = Epoch;
                    }
                    MAG_TimelyModifyMagneticModel(startdate, MagneticModels[epochs], TimedMagneticModel); /*This modifies the Magnetic coefficients to the correct date. */
                    MAG_Summation(LegendreFunction, TimedMagneticModel, *SphVariables, CoordSpherical, &MagneticResultsSph); /* Accumulate the spherical harmonic coefficients Equations 10:12 , WMM Technical report*/
                    MAG_SecVarSummation(LegendreFunction, TimedMagneticModel, *SphVariables, CoordSpherical, &MagneticResultsSphVar); /*Sum the Secular Variation Coefficients, Equations 13:15 , WMM Technical report  */
                    MAG_RotateMagneticVector(CoordSpherical, minimum, MagneticResultsSph, &MagneticResultsGeo); /* Map the computed Magnetic fields to Geodetic coordinates Equation 16 , WMM Technical report */
                    MAG_RotateMagneticVector(CoordSpherical, minimum, MagneticResultsSphVar, &MagneticResultsGeoVar); /* Map the secular variation field components to Geodetic coordinates, Equation 17 , WMM Technical report*/
                    MAG_CalculateGeoMagneticElements(&MagneticResultsGeo, &GeoMagneticElements); /* Calculate the Geomagnetic elements, Equation 18 , WMM Technical report */
                    MAG_CalculateSecularVariationElements(MagneticResultsGeoVar, &GeoMagneticElements); /*Calculate the secular variation of each of the Geomagnetic elements, Equation 19, WMM Technical report*/
                    
                    if(ElementOption>=17)
                        MAG_Gradient(Ellip, minimum, TimedMagneticModel, &Gradient);

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
                            /*16. Yearly rate of change in grid variation*/;
                            break;
                        case 17:
                            PrintElement = Gradient.GradPhi.X;
                            break;
                        case 18:
                            PrintElement = Gradient.GradPhi.Y;
                            break;
                        case 19:
                            PrintElement = Gradient.GradPhi.Z;
                            break;
                        case 20:
                            PrintElement = Gradient.GradLambda.X;
                            break;
                        case 21:
                            PrintElement = Gradient.GradLambda.Y;
                            break;
                        case 22:
                            PrintElement = Gradient.GradLambda.Z;
                            break;
                        case 23:
                            PrintElement = Gradient.GradZ.X;
                            break;
                        case 24:
                            PrintElement = Gradient.GradZ.Y;
                            break;
                        case 25:
                            PrintElement = Gradient.GradZ.Z;
                            break;
                        case 26: 
                            PrintFullGradient = 1;
                        default:
                            PrintElement = GeoMagneticElements.Decl; /* 1. Angle between the magnetic field vector and true north, positive east*/
                    }


                    if(PrintFullGradient == 1)
                    {
                        if(Geoid.UseGeoid == 1)
                        {
                            if(PrintOption == 1) fprintf(fileout, "%5.2lf %6.2lf %8.4lf %7.2lf %10.2lf %10.2lf %10.2lf %10.2lf %10.2lf %10.2lf %10.2lf %10.2lf %10.2lf\n", minimum.phi, minimum.lambda, minimum.HeightAboveGeoid, startdate.DecimalYear, 
                                    Gradient.GradPhi.X, Gradient.GradPhi.Y, Gradient.GradPhi.Z, 
                                    Gradient.GradLambda.X, Gradient.GradLambda.Y, Gradient.GradLambda.Z,
                                    Gradient.GradZ.X, Gradient.GradZ.Y, Gradient.GradZ.Z);
                            else printf("%5.2lf %6.2lf %8.4lf %7.2lf %10.2lf %10.2lf %10.2lf %10.2lf %10.2lf %10.2lf %10.2lf %10.2lf %10.2lf\n", minimum.phi, minimum.lambda, minimum.HeightAboveGeoid, startdate.DecimalYear, 
                                    Gradient.GradPhi.X, Gradient.GradPhi.Y, Gradient.GradPhi.Z, 
                                    Gradient.GradLambda.X, Gradient.GradLambda.Y, Gradient.GradLambda.Z,
                                    Gradient.GradZ.X, Gradient.GradZ.Y, Gradient.GradZ.Z);
                        } else
                        {
                            if(PrintOption == 1) fprintf(fileout, "%5.2lf %6.2lf %8.4lf %7.2lf %10.2lf\n", minimum.phi, minimum.lambda, minimum.HeightAboveEllipsoid, startdate.DecimalYear, PrintElement);
                            else printf("%5.2lf %6.2lf %8.4lf %7.2lf %10.2lf\n", minimum.phi, minimum.lambda, minimum.HeightAboveEllipsoid, startdate.DecimalYear, PrintElement);
                        }
                    } else
                    {
                    if(Geoid.UseGeoid == 1)
                    {
                        if(PrintOption == 1) fprintf(fileout, "%5.2lf %6.2lf %8.4lf %7.2lf %10.2lf\n", minimum.phi, minimum.lambda, minimum.HeightAboveGeoid, startdate.DecimalYear, PrintElement);
                        else printf("%5.2lf %6.2lf %8.4lf %7.2lf %10.2lf\n", minimum.phi, minimum.lambda, minimum.HeightAboveGeoid, startdate.DecimalYear, PrintElement);
                    } else
                    {
                        if(PrintOption == 1) fprintf(fileout, "%5.2lf %6.2lf %8.4lf %7.2lf %10.2lf\n", minimum.phi, minimum.lambda, minimum.HeightAboveEllipsoid, startdate.DecimalYear, PrintElement);
                        else printf("%5.2lf %6.2lf %8.4lf %7.2lf %10.2lf\n", minimum.phi, minimum.lambda, minimum.HeightAboveEllipsoid, startdate.DecimalYear, PrintElement);
                    }
                    }
                        /*if(PrintOption == 1) fprintf(fileout, "%6.2lf %6.2lf %10.2lf\n", minimum.lambda, minimum.phi, PrintElement);
                        else printf("%6.2lf %6.2lf %10.2lf\n", minimum.lambda, minimum.phi, PrintElement);*/

                } /* year loop */

            } /*Longitude Loop */

        } /* Latitude Loop */

    } /* Altitude Loop */
    if(PrintOption == 1) fclose(fileout);

    for(i = 0; i < epochs + 1; i++) MAG_FreeMagneticModelMemory(MagneticModels[i]);
    MAG_FreeMagneticModelMemory(TimedMagneticModel);
    MAG_FreeLegendreMemory(LegendreFunction);
    MAG_FreeSphVarMemory(SphVariables);



    printf("\nPress any key to exit...\n");
    getchar();

    return 0;
}

