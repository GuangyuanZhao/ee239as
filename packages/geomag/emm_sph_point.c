//---------------------------------------------------------------------------

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>


#include "GeomagnetismHeader.h"
//#include "UnreleasedFunctions.h"
#include "EGM9615.h"
//#include "GeomagnetismLibrary.c"



//---------------------------------------------------------------------------

/* The Geomagnetism Library is used to make a command prompt program. The program prompts
the user to enter a location, performs the computations and prints the results to the
standard output. The program expects the files GeomagnetismLibrary.c, GeomagnetismHeader.h,
EMM2015.COF, EMM2015SV.COF and EGM9615.h to be in the same directory. 

Manoj.C.Nair
Nov 23, 2009

 *  Revision Number: $Revision: 1324 $
 *  Last changed by: $Author: awoods $
 *  Last changed on: $Date: 2015-05-13 16:04:23 -0600 (Wed, 13 May 2015) $

 */

int main()
{

    MAGtype_Ellipsoid Ellip;
    MAGtype_CoordSpherical CoordSpherical;
    MAGtype_CoordGeodetic CoordGeodetic;
    MAGtype_Date UserDate;
    MAGtype_GeoMagneticElements GeoMagneticElements;
    MAGtype_Gradient Gradient;
    MAGtype_Geoid Geoid;
    char ans[20], b;
    char filename[] = "EMM2015.COF";
    char filenameSV[] = "EMM2015SV.COF";
    int NumTerms, Epoch, LoadedEpoch = -1, i, Flag = 1, nMax = 0, nMaxEMM, index;
    const int epochs = 16;
    char VersionDate_Large[] = "$Date: 2015-05-13 16:04:23 -0600 (Wed, 13 May 2015) $";
    char VersionDate[12];
    /*The main field for the Magnetic Models 2000-2015 are stored in MagneticModel[0] to
     MagneticModel[15].  The full model including the crustal field is in MagneticModel[16]*/
    MAGtype_MagneticModel *MagneticModels[17], *TimedMagneticModel;
    /* Memory allocation */

    Gradient.UseGradient = 1; 
    strncpy(VersionDate, VersionDate_Large + 39, 11);
    VersionDate[11] = '\0';
    for(Epoch = epochs; Epoch >= 0; Epoch--)
    {
        if(Epoch == epochs-1)
            Epoch--;
        if(Epoch < epochs-1)
        {
            sprintf(filename, "EMM%d.COF", Epoch + 2000);
            sprintf(filenameSV, "EMM%dSV.COF", Epoch + 2000);
        }

        if(!MAG_robustReadMagneticModel_Large(filename, filenameSV, &MagneticModels[Epoch]))  {
            printf("\n EMM%d.COF or EMM%dSV.COF not found.  Press enter to exit... \n ", Epoch+2000, Epoch+2000);
            fgets(ans, 20, stdin);
            return 1;
        }
        MagneticModels[Epoch]->CoefficientFileEndDate = MagneticModels[epochs]->CoefficientFileEndDate;
    }
    /*Create Main Field Model for EMM2015*/
    nMaxEMM = MagneticModels[0]->nMax;
    NumTerms = ((nMaxEMM + 1) * (nMaxEMM + 2) / 2);
    MagneticModels[epochs - 1] = MAG_AllocateModelMemory(NumTerms);
    MagneticModels[epochs - 1]->nMax = nMaxEMM;
    MagneticModels[epochs - 1]->nMaxSecVar = nMaxEMM;
    MagneticModels[epochs - 1]->epoch = MagneticModels[0]->epoch + epochs - 1;
    MAG_AssignMagneticModelCoeffs(MagneticModels[epochs - 1], MagneticModels[epochs], MagneticModels[epochs - 1]->nMax, MagneticModels[epochs - 1]->nMaxSecVar);


    /*Allocate Memory for TimedMagneticModel*/

    nMax = MagneticModels[epochs]->nMax;
    NumTerms = ((nMax + 1) * (nMax + 2) / 2);
    TimedMagneticModel = MAG_AllocateModelMemory(NumTerms); /* For storing the time modified WMM Model parameters */
    for(i = 0; i < epochs; i++) if(MagneticModels[i] == NULL || TimedMagneticModel == NULL)
        {
            MAG_Error(2);
        }
    MAG_SetDefaults(&Ellip, &Geoid); /* Set default values and constants */
    /* Check for Geographic Poles */


    //MAG_InitializeGeoid(&Geoid);    /* Read the Geoid file */
    /* Set EGM96 Geoid parameters */
    Geoid.GeoidHeightBuffer = GeoidHeightBuffer;
    Geoid.Geoid_Initialized = 1;
    /* Set EGM96 Geoid parameters END */
    b = MAG_GeomagIntroduction_EMM(MagneticModels[epochs-1], VersionDate);

    while(Flag == 1 && b != 'x')
    {
        if(MAG_GetUserInput(MagneticModels[0], &Geoid, &CoordGeodetic, &UserDate) == 1) /*Get User Input */
        {
            Epoch = ((int) UserDate.DecimalYear - MagneticModels[0]->epoch);
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
            MAG_GeodeticToSpherical(Ellip, CoordGeodetic, &CoordSpherical); /*Convert from geodetic to Spherical Equations: 17-18, WMM Technical report*/
            MAG_TimelyModifyMagneticModel(UserDate, MagneticModels[epochs], TimedMagneticModel); /* Time adjust the coefficients, Equation 19, WMM Technical report */
            MAG_Geomag(Ellip, CoordSpherical, CoordGeodetic, TimedMagneticModel, &GeoMagneticElements); /* Computes the geoMagnetic field elements and their time change*/
            MAG_CalculateGridVariation(CoordGeodetic, &GeoMagneticElements);
            MAG_PrintUserData(GeoMagneticElements, CoordGeodetic, UserDate, TimedMagneticModel, &Geoid); /* Print the results */
            printf("\nWould you like to print the gradient? (y or n)\n");
            fgets(ans, 20, stdin);
            switch(ans[0]) {
                case 'Y':
                case 'y':
                    MAG_Gradient(Ellip, CoordGeodetic, TimedMagneticModel, &Gradient);
                    MAG_PrintGradient(Gradient);
                    break;
                case 'N':
                case 'n':
                default:
                    break;
            }
        }
        printf("\n\n Do you need more point data ? (y or n) \n ");
        strcpy(ans,"");
        fgets(ans, 20, stdin);
        switch(ans[0]) {
            case 'Y':
            case 'y':
                Flag = 1;
                break;
            case 'N':
            case 'n':
                Flag = 0;
                break;
            default:
                Flag = 0;
                break;
        }

    }


    MAG_FreeMagneticModelMemory(TimedMagneticModel);
    for(i = 0; i < epochs + 1; i++) MAG_FreeMagneticModelMemory(MagneticModels[i]);



    return 0;
}
