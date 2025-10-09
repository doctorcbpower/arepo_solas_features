#include "config.h"
#include "SNIIYields.h"
#include "../data/Portinari+1998/Portinari+1998_Struct.h"
#include "../data/Nomoto+2013/Nomoto+2013_Struct.h"
#include "../data/Nomoto+2013/Nomoto+2013HN_Struct.h"

/*! \file SNIIYields.c
 * \brief This file has functions regarding SNII yields.
 */

/*!
 * Following variables hold the grid size of the adopted yields table and values
 * of metallicity.
 */
int CELibSNIIYields_Metallicity; 
int CELibSNIIYields_Mass; 
double *CELibSNIIYieldsZ; 

struct CELibStructSNIIYields *pCELibSNIIYields; 
#define CELibSNIIYields(x,y) (pCELibSNIIYields[CELibSNIIYields_Mass*(x)+(y)])

/*!
 * Array of structure for Nomoto+2013 yields. Since the yields table includes
 * Pop III yields, several special cares are necessary. For this, this array of
 * structure is used.
 */
struct CELibStructSNIIYields *pCELibSNIIYieldsN13 = NULL; 

/*
 * This is the IMF weighted, integrated yield tables.
 */
struct CELibStructSNIIYields *CELibSNIIYieldsIntegrated = NULL;


static void pCELibDumpSNIIYieldsInterpolatedValues(const char OutDir[]);
static void pCElibMakeSNIIYieldsIntegratedTable(void);
static void pCElibWriteSNIIYieldsIntegratedTable(const char OutDir[]);
static void pCElibWriteSNIIYieldsReturnMassFraction(const char OutDir[]);
static void pCElibWriteSNIIYieldsHyperNovaYields(const char OutDir[]);
static void pCELibCopySNIIYieldsRanges(void);


static char SNIIYieldsElementName[CELibYield_Number][MaxCharactersInLine]; 

/*!
 * Set yield names.
 */
static void pCELibSNIIYieldsSetElementName(void){

    strcpy(SNIIYieldsElementName[CELibYield_H],    "H"); 
    strcpy(SNIIYieldsElementName[CELibYield_He],   "He"); 
    strcpy(SNIIYieldsElementName[CELibYield_C],    "C"); 
    strcpy(SNIIYieldsElementName[CELibYield_N],    "N"); 
    strcpy(SNIIYieldsElementName[CELibYield_O],    "O"); 
    strcpy(SNIIYieldsElementName[CELibYield_Ne],   "Ne"); 
    strcpy(SNIIYieldsElementName[CELibYield_Mg],   "Mg"); 
    strcpy(SNIIYieldsElementName[CELibYield_Si],   "Si"); 
    strcpy(SNIIYieldsElementName[CELibYield_S],    "S"); 
    strcpy(SNIIYieldsElementName[CELibYield_Ca],   "Ca"); 
    strcpy(SNIIYieldsElementName[CELibYield_Fe],   "Fe"); 
    strcpy(SNIIYieldsElementName[CELibYield_Ni],   "Ni"); 
    strcpy(SNIIYieldsElementName[CELibYield_Eu],   "Eu"); 

    return ;
}


/*!
 * Dump all data.
 */
static void pCELibDumpSNIIYieldsAllData(char OutDir[]){

    if(CELibRunParameters.TestMode == false)
        return ;

    FILE *fp;
    char fname[MaxCharactersInLine];
    if(CELibRunParameters.SNIIYieldsTableID == CELibSNIIYieldsTableID_P98){
        Snprintf(fname,"%s/DumpYieldsP98.dat",OutDir);
    } else if(CELibRunParameters.SNIIYieldsTableID == CELibSNIIYieldsTableID_N13){
        Snprintf(fname,"%s/DumpYieldsN13.dat",OutDir);
    }
    FileOpen(fp,fname,"w");
    for(int i=0;i<CELibSNIIYields_Metallicity;i++){
        for(int k=0;k<CELibSNIIYields_Mass;k++){
            fprintf(fp,"###### (%d,%d)\n",i,k);
            fprintf(fp,"Mass %g\n",CELibSNIIYields(i,k).Mass);
            fprintf(fp,"Mrem %g\n",CELibSNIIYields(i,k).Mr);
            fprintf(fp,"Msw %g\n",CELibSNIIYields(i,k).Msw);
            fprintf(fp,"Msn %g\n",CELibSNIIYields(i,k).Msn);
            fprintf(fp,"Energy %g\n",CELibSNIIYields(i,k).Energy);
        }
    }
    fclose(fp);
    
    return;
}

/*!
 * Dump all integrated data.
 */
static void pCELibDumpSNIIYieldsAllIntegratedData(char OutDir[]){

    if(CELibRunParameters.TestMode == false)
        return ;

    FILE *fp;
    char fname[MaxCharactersInLine];
    if(CELibRunParameters.SNIIYieldsTableID == CELibSNIIYieldsTableID_P98){
        Snprintf(fname,"%s/DumpYieldsIntegratedP98.dat",OutDir);
    } else if(CELibRunParameters.SNIIYieldsTableID == CELibSNIIYieldsTableID_N13){
        Snprintf(fname,"%s/DumpYieldsIntegratedN13.dat",OutDir);
    }
    FileOpen(fp,fname,"w");
    for(int i=0;i<CELibSNIIYields_Metallicity;i++){
        fprintf(fp,"###### %d\n",i);
        fprintf(fp,"Mass %g\n",CELibSNIIYieldsIntegrated[i].Mass);
        fprintf(fp,"Mrem %g\n",CELibSNIIYieldsIntegrated[i].Mr);
        fprintf(fp,"Msw %g\n",CELibSNIIYieldsIntegrated[i].Msw);
        fprintf(fp,"Msn %g\n",CELibSNIIYieldsIntegrated[i].Msn);
        fprintf(fp,"Energy %g\n",CELibSNIIYieldsIntegrated[i].Energy);
    }
    fclose(fp);
    
    return;
}

/*!
 * This function returns the minimum stellar mass which becomes a type II SN.
 */
double CELibGetSNIIYieldsMinExplosionMass(const int IndexMetal){

    for(int k=0;k<CELibSNIIYields_Mass;k++){
        if(CELibSNIIYields(IndexMetal,k).Energy > 0.e0){
            return CELibSNIIYields(IndexMetal,k).Mass;
        }
    }

    return -1;
}

/*!
 * This function returns the maximum stellar mass which becomes a type II SN.
 */
double CELibGetSNIIYieldsMaxExplosionMass(const int IndexMetal){

    for(int k=0;k<CELibSNIIYields_Mass;k++){
        if(CELibSNIIYields(IndexMetal,CELibSNIIYields_Mass-1-k).Energy > 0.e0){
            return CELibSNIIYields(IndexMetal,CELibSNIIYields_Mass-1-k).Mass;
        }
    }

    return -1;
}


/*!
 * This function shows both the minimum and maximum stellar masses which become
 * type II SNe. 
 */
static void pCELibShowSNIIYieldsExplosionMassRange(void){

    fprintf(stderr,"//\t Show SNII explosion mass range.\n");
    for(int i=0;i<CELibSNIIYields_Metallicity;i++){
        double Max = CELibGetSNIIYieldsMaxExplosionMass(i);
        double Min = CELibGetSNIIYieldsMinExplosionMass(i);
        fprintf(stderr,"//\t\t Z= %g,  %g -- %g\n",
                CELibSNIIYields(i,0).Metallicity,Min,Max);
    }
    return ;
}


/*!
 * This function retuns hyprenova yields, the remenant mass, the ejecta mass via
 * a SN and stellar winds, and the released energy.
 */
static int pCELibGetSNIIYieldsHyperNovaYields(const int MetalID, const double Mass, double Elements[], double *Mr, double *Msw, double *Msn, double *Energy){

    int IndexMax;
    if(MetalID == 0){
        IndexMax = CELibSNIIYields_Mass_N13HN_Zero;
    } else {
        IndexMax = CELibSNIIYields_Mass_N13HN_Normal;
    }

    if((Mass < CELibSNIIYieldsN13HN[MetalID][0].Mass)||(CELibSNIIYieldsN13HN[MetalID][IndexMax-1].Mass<Mass)){
        memset(Elements,0,sizeof(double)*CELibYield_Number);
        return 0;
    }

    int IndexMass = 0;
    for(int i=1;i<IndexMax;i++){
        if(CELibSNIIYieldsN13HN[MetalID][i].Mass > Mass){
            IndexMass = i-1;
            break;
        }
    }

    double MassMin = CELibSNIIYieldsN13HN[MetalID][IndexMass].Mass;
    double MassMax = CELibSNIIYieldsN13HN[MetalID][IndexMass+1].Mass;
    double InvdMass = 1.0/(MassMax-MassMin);
    for(int i=0;i<CELibYield_Number;i++){
        Elements[i] = 
            (CELibSNIIYieldsN13HN[MetalID][IndexMass+1].Elements[i]-CELibSNIIYieldsN13HN[MetalID][IndexMass].Elements[i])*InvdMass
            *(Mass-CELibSNIIYieldsN13HN[MetalID][IndexMass].Mass)+CELibSNIIYieldsN13HN[MetalID][IndexMass].Elements[i];
    }

    *Mr = (CELibSNIIYieldsN13HN[MetalID][IndexMass+1].Mr-CELibSNIIYieldsN13HN[MetalID][IndexMass].Mr)*InvdMass
        *(Mass-CELibSNIIYieldsN13HN[MetalID][IndexMass].Mass)+CELibSNIIYieldsN13HN[MetalID][IndexMass].Mr;
    *Msn = (CELibSNIIYieldsN13HN[MetalID][IndexMass+1].Msn-CELibSNIIYieldsN13HN[MetalID][IndexMass].Msn)*InvdMass
        *(Mass-CELibSNIIYieldsN13HN[MetalID][IndexMass].Mass)+CELibSNIIYieldsN13HN[MetalID][IndexMass].Msn;
    *Msw = Mass - *Mr -*Msn;
    *Energy = (CELibSNIIYieldsN13HN[MetalID][IndexMass+1].Energy-CELibSNIIYieldsN13HN[MetalID][IndexMass].Energy)*InvdMass
        *(Mass-CELibSNIIYieldsN13HN[MetalID][IndexMass].Mass)+CELibSNIIYieldsN13HN[MetalID][IndexMass].Energy;

    return 1;
}

/*!
 * This function blends hypernova yields with the oridinary type II SNe yields.
 * The hyprenova faction is defined in CELibRunParameters.SNIIHyperNovaFraction.
 */
static void pCELibBlendSNIIYieldsHN(struct CELibStructSNIIYields *y, const int MetalID){

    double HNElements[CELibYield_Number]; // Blend hypernova
    double Mr, Msw, Msn, Energy;
    if(pCELibGetSNIIYieldsHyperNovaYields(MetalID,y->Mass,HNElements,&Mr,&Msw,&Msn,&Energy)){
        y->Mr *= (1.0-CELibRunParameters.SNIIHyperNovaFraction);
        y->Mr += CELibRunParameters.SNIIHyperNovaFraction*Mr;
        y->Msn *= (1.0-CELibRunParameters.SNIIHyperNovaFraction);
        y->Msn += CELibRunParameters.SNIIHyperNovaFraction*Msn;
        y->Msw *= (1.0-CELibRunParameters.SNIIHyperNovaFraction);
        y->Msw += CELibRunParameters.SNIIHyperNovaFraction*Msw;
        y->Energy *= (1.0-CELibRunParameters.SNIIHyperNovaFraction);
        y->Energy += CELibRunParameters.SNIIHyperNovaFraction*Energy;

        for(int l=0;l<CELibYield_Number;l++){
            y->Elements[l] *= (1.0-CELibRunParameters.SNIIHyperNovaFraction);
            y->Elements[l] += CELibRunParameters.SNIIHyperNovaFraction*HNElements[l];
        }
    }

    return ;
}

/*!
 * This function initializes the yield table of type II SNe. 
 */
void CELibInitSNIIYields(void){

    if(CELibRunParameters.SNIIYieldsTableID == CELibSNIIYieldsTableID_P98){
        CELibSNIIYields_Metallicity = CELibSNIIYields_Metallicity_P98; 
        CELibSNIIYields_Mass = CELibSNIIYields_Mass_P98;  
        pCELibSNIIYields = CELibSNIIYieldsP98[0]; 
        CELibSNIIYieldsIntegrated = realloc(CELibSNIIYieldsIntegrated,
                sizeof(struct CELibStructSNIIYields)*CELibSNIIYields_Metallicity);

    } else if(CELibRunParameters.SNIIYieldsTableID == CELibSNIIYieldsTableID_N13){
        CELibSNIIYields_Mass = CELibSNIIYields_Mass_N13; 
        CELibSNIIYields_Metallicity = CELibSNIIYields_Metallicity_N13; 
        
        if(CELibRunParameters.PopIIISNe == 1){ // Pop III mode on
            CELibSNIIYields_Metallicity ++;

            pCELibSNIIYieldsN13 = realloc(pCELibSNIIYieldsN13,
                    sizeof(struct CELibStructSNIIYields)*CELibSNIIYields_Metallicity*CELibSNIIYields_Mass);
            for(int k=0;k<CELibSNIIYields_Mass;k++){
                pCELibSNIIYieldsN13[k] = CELibSNIIYieldsN13[0][k];
                pCELibSNIIYieldsN13[k].Metallicity = CELibRunParameters.PopIIIMetallicity;
                pCELibSNIIYieldsN13[CELibSNIIYields_Mass+k] = CELibSNIIYieldsN13[1][k];
                pCELibSNIIYieldsN13[CELibSNIIYields_Mass+k].Metallicity = 1.1*CELibRunParameters.PopIIIMetallicity;


                pCELibBlendSNIIYieldsHN(pCELibSNIIYieldsN13+k,0);
                pCELibBlendSNIIYieldsHN(pCELibSNIIYieldsN13+CELibSNIIYields_Mass+k,0);
            }
            for(int i=1;i<CELibSNIIYields_Metallicity_N13;i++){
                for(int k=0;k<CELibSNIIYields_Mass;k++){
                    pCELibSNIIYieldsN13[(i+1)*CELibSNIIYields_Mass+k] = CELibSNIIYieldsN13[i][k];

                    pCELibBlendSNIIYieldsHN(pCELibSNIIYieldsN13+(i+1)*CELibSNIIYields_Mass+k,i);
                }
            }
        } else {  // Pop III mode off
            pCELibSNIIYieldsN13 = realloc(pCELibSNIIYieldsN13,
                    sizeof(struct CELibStructSNIIYields)*CELibSNIIYields_Metallicity*CELibSNIIYields_Mass);
            for(int k=0;k<CELibSNIIYields_Mass;k++){
                pCELibSNIIYieldsN13[k] = CELibSNIIYieldsN13[1][k];

                pCELibBlendSNIIYieldsHN(pCELibSNIIYieldsN13+k,0);
            }
            for(int i=1;i<CELibSNIIYields_Metallicity;i++){
                for(int k=0;k<CELibSNIIYields_Mass;k++){
                    pCELibSNIIYieldsN13[i*CELibSNIIYields_Mass+k] = CELibSNIIYieldsN13[i][k];
                    pCELibBlendSNIIYieldsHN(pCELibSNIIYieldsN13+i*CELibSNIIYields_Mass+k,i);
                }
            }
        }

        pCELibSNIIYields = pCELibSNIIYieldsN13; 

        CELibSNIIYieldsIntegrated = realloc(CELibSNIIYieldsIntegrated,
                sizeof(struct CELibStructSNIIYields)*CELibSNIIYields_Metallicity);
    }

    CELibSNIIYieldsZ = realloc(CELibSNIIYieldsZ,sizeof(double)*CELibSNIIYields_Metallicity);
    for(int i=0;i<CELibSNIIYields_Metallicity;i++){
        CELibSNIIYieldsZ[i] = CELibSNIIYields(i,0).Metallicity;
    }

    pCELibSNIIYieldsSetElementName();
    pCElibMakeSNIIYieldsIntegratedTable();

    if(CELibRunParameters.TestMode){
        pCElibWriteSNIIYieldsIntegratedTable("./CELib/CELibSNII");
        pCELibDumpSNIIYieldsInterpolatedValues("./CELib/CELibSNII");

        pCELibDumpSNIIYieldsAllData("./CELib/CELibSNII");
        pCELibDumpSNIIYieldsAllIntegratedData("./CELib/CELibSNII");
        pCElibWriteSNIIYieldsReturnMassFraction("./CELib/CELibSNII");
        pCElibWriteSNIIYieldsHyperNovaYields("./CELib/CELibSNII");
        pCELibShowSNIIYieldsExplosionMassRange();
    }

    pCELibCopySNIIYieldsRanges();

    return ;
}

/*! 
 * This function retuns an interpolated yeild of the element "TargetElement"
 * from a star of which mass is "Mass".
 */
static double pCELibGetSNIIYieldsInterpolatedYieldMass(const int IndexMetal, const int TargetElement, const double Mass){

    if(CELibSNIIYields(IndexMetal,0).Mass > Mass)
        return 0.e0;
    if(CELibSNIIYields(IndexMetal,CELibSNIIYields_Mass-1).Mass < Mass)
        return 0.e0;

    int IndexMass = 0;
    for(int i=1;i<CELibSNIIYields_Mass;i++){
        if(CELibSNIIYields(IndexMetal,i).Mass > Mass){
            IndexMass = i-1;
            break;
        }
    }
    
    return (CELibSNIIYields(IndexMetal,IndexMass+1).Elements[TargetElement]-CELibSNIIYields(IndexMetal,IndexMass).Elements[TargetElement])/
        (CELibSNIIYields(IndexMetal,IndexMass+1).Mass-CELibSNIIYields(IndexMetal,IndexMass).Mass)*
        (Mass-CELibSNIIYields(IndexMetal,IndexMass).Mass)+CELibSNIIYields(IndexMetal,IndexMass).Elements[TargetElement];
}

/*!
 *  This function retuns an interpolated remnant mass from a given mass star.
 */
static double pCELibGetSNIIYieldsInterpolatedRemnantMass(const int IndexMetal, const double Mass){

    int IndexMass = 0;
    if(CELibSNIIYields(IndexMetal,0).Mass > Mass)
        return 0.e0;
    if(CELibSNIIYields(IndexMetal,CELibSNIIYields_Mass-1).Mass < Mass)
        return 0.e0;

    for(int i=1;i<CELibSNIIYields_Mass;i++){
        if(CELibSNIIYields(IndexMetal,i).Mass > Mass){
            IndexMass = i-1;
            break;
        }
    }
    
    return (CELibSNIIYields(IndexMetal,IndexMass+1).Mr-CELibSNIIYields(IndexMetal,IndexMass).Mr)/
        (CELibSNIIYields(IndexMetal,IndexMass+1).Mass-CELibSNIIYields(IndexMetal,IndexMass).Mass)*
        (Mass-CELibSNIIYields(IndexMetal,IndexMass).Mass)+CELibSNIIYields(IndexMetal,IndexMass).Mr;
}

/*!
 * This function retuns an interpolated mass released via stellar winds from a
 * given mass star.
 */
static double pCELibGetSNIIYieldsInterpolatedStellarWindMass(const int IndexMetal, const double Mass){

    int IndexMass = 0;
    if(CELibSNIIYields(IndexMetal,0).Mass > Mass)
        return 0.e0;
    if(CELibSNIIYields(IndexMetal,CELibSNIIYields_Mass-1).Mass < Mass)
        return 0.e0;

    for(int i=1;i<CELibSNIIYields_Mass;i++){
        if(CELibSNIIYields(IndexMetal,i).Mass > Mass){
            IndexMass = i-1;
            break;
        }
    }
    
    return (CELibSNIIYields(IndexMetal,IndexMass+1).Msw-CELibSNIIYields(IndexMetal,IndexMass).Msw)/
        (CELibSNIIYields(IndexMetal,IndexMass+1).Mass-CELibSNIIYields(IndexMetal,IndexMass).Mass)*
        (Mass-CELibSNIIYields(IndexMetal,IndexMass).Mass)+CELibSNIIYields(IndexMetal,IndexMass).Msw;
}

/*! 
 * This function retuns an interpolated mass released via SN from a given mass
 * star.
 */
static double pCELibGetSNIIYieldsInterpolatedSuperNovaMass(const int IndexMetal, const double Mass){

    int IndexMass = 0;
    if(CELibSNIIYields(IndexMetal,0).Mass > Mass)
        return 0.e0;
    if(CELibSNIIYields(IndexMetal,CELibSNIIYields_Mass-1).Mass < Mass)
        return 0.e0;

    for(int i=1;i<CELibSNIIYields_Mass;i++){
        if(CELibSNIIYields(IndexMetal,i).Mass > Mass){
            IndexMass = i-1;
            break;
        }
    }
    
    return (CELibSNIIYields(IndexMetal,IndexMass+1).Msn-CELibSNIIYields(IndexMetal,IndexMass).Msn)/
        (CELibSNIIYields(IndexMetal,IndexMass+1).Mass-CELibSNIIYields(IndexMetal,IndexMass).Mass)*
        (Mass-CELibSNIIYields(IndexMetal,IndexMass).Mass)+CELibSNIIYields(IndexMetal,IndexMass).Msn;
}

/*! 
 * This function retuns the interpolated energy released from a given mass star.
 */
static double pCELibGetSNIIYieldsInterpolatedEnergy(const int IndexMetal, const double Mass){

    int IndexMass = 0;
    if(CELibSNIIYields(IndexMetal,0).Mass > Mass)
        return 0.e0;
    if(CELibSNIIYields(IndexMetal,CELibSNIIYields_Mass-1).Mass < Mass)
        return 0.e0;

    for(int i=1;i<CELibSNIIYields_Mass;i++){
        if(CELibSNIIYields(IndexMetal,i).Mass > Mass){
            IndexMass = i-1;
            break;
        }
    }
    
    return (CELibSNIIYields(IndexMetal,IndexMass+1).Energy-CELibSNIIYields(IndexMetal,IndexMass).Energy)/
        (CELibSNIIYields(IndexMetal,IndexMass+1).Mass-CELibSNIIYields(IndexMetal,IndexMass).Mass)*
        (Mass-CELibSNIIYields(IndexMetal,IndexMass).Mass)+CELibSNIIYields(IndexMetal,IndexMass).Energy;
}

/*!
 * This function writes the net yields (y) and the stellar yields (Y) as a
 * function of the progenitor stellar mass in a file. For the stellar yields,
 * projenitors' abundance pattern is the same as the solar one. The linear
 * interpolation is employed.
 */
static void pCELibDumpSNIIYieldsInterpolatedValues(const char OutDir[]){

    const int NSample = 200;
    MakeDir(OutDir);

    for(int i=0;i<CELibSNIIYields_Metallicity;i++){
        for(int j=0;j<CELibYield_Number;j++){
            FILE *fp;
            char fname[MaxCharactersInLine];
            Snprintf(fname,"%s/SNIIYieldsInterpolated.%02d.%02d",OutDir,i,j);
            FileOpen(fp,fname,"w");

            fprintf(fp,"#Metal ID = %d, TargetID = %d %s\n",i,j,SNIIYieldsElementName[j]);
            fprintf(fp,"#Metallicity = %g\n",CELibSNIIYields(i,0).Metallicity);
            double dMass = (CELibSNIIYields(i,CELibSNIIYields_Mass-1).Mass-CELibSNIIYields(i,0).Mass)/NSample;
            for(int k=0;k<NSample;k++){
                double Mass = dMass*k + CELibSNIIYields(i,0).Mass;
                double Value = pCELibGetSNIIYieldsInterpolatedYieldMass(i,j,Mass);
                double FeValue = pCELibGetSNIIYieldsInterpolatedYieldMass(i,CELibYield_Fe,Mass);
                fprintf(fp,"%g %g %g %g\n",Mass,Value,Value/FeValue,log10(Value/FeValue));
            }
            fclose(fp);
        }
        {
            FILE *fp;
            char fname[MaxCharactersInLine];
            // Mr, Msw and Energy interporlated values.
            Snprintf(fname,"%s/SNIIMrMswEInterpolated.%02d",OutDir,i);
            FileOpen(fp,fname,"w");

            fprintf(fp,"#Metal ID = %d\n",i);
            fprintf(fp,"#Metallicity = %g\n",CELibSNIIYields(i,0).Metallicity);
            fprintf(fp,"#Mass Mr Msw Msn Mejecta E\n");
            double dMass = (CELibSNIIYields(i,CELibSNIIYields_Mass-1).Mass-CELibSNIIYields(i,0).Mass)/NSample;
            for(int k=0;k<NSample;k++){
                double Mass = dMass*k + CELibSNIIYields(i,0).Mass;
                fprintf(fp,"%g %g %g %g %g %g\n",Mass,pCELibGetSNIIYieldsInterpolatedRemnantMass(i,Mass),
                        pCELibGetSNIIYieldsInterpolatedStellarWindMass(i,Mass),
                        pCELibGetSNIIYieldsInterpolatedSuperNovaMass(i,Mass),
                        pCELibGetSNIIYieldsInterpolatedStellarWindMass(i,Mass)+
                        pCELibGetSNIIYieldsInterpolatedSuperNovaMass(i,Mass),
                        pCELibGetSNIIYieldsInterpolatedEnergy(i,Mass));
            }
            fclose(fp);
        }
    }


    return ;
}

static int pCELibSNIIYieldsCurrentIMFType;
static int pCELibSNIIYieldsCurrentMetallicity;
static int pCELibSNIIYieldsCurrentElement;

/*!
 * IMF weighted remnant mass.
 */
static double pCELibGetSNIIYieldsMrFunction(const double Mass){
    double Mr = pCELibGetSNIIYieldsInterpolatedRemnantMass(pCELibSNIIYieldsCurrentMetallicity,Mass);
    return Mr*CELibIMF[pCELibSNIIYieldsCurrentIMFType].IMFFunctionPerMass(Mass);
}

/*!
 * This function returns an IMF weighted remnant mass within a given mass range.
 */
static double pCELibGetSNIIYieldsMrInGivenMassRange(const double LowerMass_in_SolarMass, const double UpperMass_in_SolarMass){
    return IntegralSimpson(LowerMass_in_SolarMass,UpperMass_in_SolarMass,CELibRunParameters.IntegrationSteps,&pCELibGetSNIIYieldsMrFunction);
}

/*!
 *
 */
static double pCELibGetSNIIYieldsYieldMassFunction(const double Mass){
    double ElementMass = pCELibGetSNIIYieldsInterpolatedYieldMass(pCELibSNIIYieldsCurrentMetallicity,pCELibSNIIYieldsCurrentElement,Mass);
    return ElementMass*CELibIMF[pCELibSNIIYieldsCurrentIMFType].IMFFunctionPerMass(Mass);
}

/*!
 *
 */
static double pCELibGetSNIIYieldsYieldMassInGivenMassRange(const double LowerMass_in_SolarMass, const double UpperMass_in_SolarMass){
    return IntegralSimpson(LowerMass_in_SolarMass,UpperMass_in_SolarMass,CELibRunParameters.IntegrationSteps,&pCELibGetSNIIYieldsYieldMassFunction);
}

/*!
 * IMF weighted mass released via stellar winds.
 */
static double pCELibGetSNIIYieldsMswFunction(const double Mass){
    double Msw = pCELibGetSNIIYieldsInterpolatedStellarWindMass(pCELibSNIIYieldsCurrentMetallicity,Mass);
    return Msw*CELibIMF[pCELibSNIIYieldsCurrentIMFType].IMFFunctionPerMass(Mass);
}

/*!
 * This function returns an IMF weighted mass released via stellar winds within a given mass range.
 */
static double pCELibGetSNIIYieldsMswInGivenMassRange(const double LowerMass_in_SolarMass, const double UpperMass_in_SolarMass){
    return IntegralSimpson(LowerMass_in_SolarMass,UpperMass_in_SolarMass,CELibRunParameters.IntegrationSteps,&pCELibGetSNIIYieldsMswFunction);
}

/*!
 * IMF weighted mass released via supernova.
 */
static double pCELibGetSNIIYieldsMsnFunction(const double Mass){
    double Msw = pCELibGetSNIIYieldsInterpolatedSuperNovaMass(pCELibSNIIYieldsCurrentMetallicity,Mass);
    return Msw*CELibIMF[pCELibSNIIYieldsCurrentIMFType].IMFFunctionPerMass(Mass);
}

/*!
 * This function returns an IMF weighted mass released via supernovae within a given mass range.
 */
static double pCELibGetSNIIYieldsMsnInGivenMassRange(const double LowerMass_in_SolarMass, const double UpperMass_in_SolarMass){
    return IntegralSimpson(LowerMass_in_SolarMass,UpperMass_in_SolarMass,CELibRunParameters.IntegrationSteps,&pCELibGetSNIIYieldsMsnFunction);
}

/*!
 * IMF weighted energy.
 */
static double pCELibGetSNIIYieldsEnergyFunction(const double Mass){
    double Energy = pCELibGetSNIIYieldsInterpolatedEnergy(pCELibSNIIYieldsCurrentMetallicity,Mass);
    return Energy*CELibIMF[pCELibSNIIYieldsCurrentIMFType].IMFFunctionPerMass(Mass);
}

/*!
 * This function returns an IMF weighted released energy within a given mass range.
 */
static double pCELibGetSNIIYieldsEnergyInGivenMassRange(const double LowerMass_in_SolarMass, const double UpperMass_in_SolarMass){
    return IntegralSimpson(LowerMass_in_SolarMass,UpperMass_in_SolarMass,CELibRunParameters.IntegrationSteps,&pCELibGetSNIIYieldsEnergyFunction);
}

/*!
 * This function makes the integrated yields data in files.
 */
static void pCElibMakeSNIIYieldsIntegratedTable(void){

    pCELibSNIIYieldsCurrentIMFType = CELibRunParameters.IMFType;

    if(CELibRunParameters.TestMode){
        if(CELibRunParameters.SNIIYieldsTableID == CELibSNIIYieldsTableID_P98){
            fprintf(stderr,"//\t CELib integrated SNII Yields Table (Portinari+1998).\n");
        } else if(CELibRunParameters.SNIIYieldsTableID == CELibSNIIYieldsTableID_N13){
            fprintf(stderr,"//\t CELib integrated SNII Yields Table (Nomoto+2013).\n");
        }
    }

    for(int i=0;i<CELibSNIIYields_Metallicity;i++){
        double MassL = CELibRunParameters.SNIILowerMass;
        double MassU = CELibRunParameters.SNIIUpperMass;

        if(CELibRunParameters.SNIIYieldsTableID == CELibSNIIYieldsTableID_N13){
            if(CELibRunParameters.PopIIISNe == 1){
                if(CELibSNIIYields(i,0).Metallicity <= CELibRunParameters.PopIIIMetallicity){
                    if(CELibRunParameters.PopIIIIMF == 1){
                        pCELibSNIIYieldsCurrentIMFType = CELibIMF_Susa;
                        MassU = CELibIMF[CELibIMF_Susa].MassMax;
                    } else {
                        pCELibSNIIYieldsCurrentIMFType = CELibRunParameters.IMFType;
                        MassU = CELibIMF[CELibRunParameters.IMFType].MassMax;
                    }
                } else {
                    pCELibSNIIYieldsCurrentIMFType = CELibRunParameters.IMFType;
                }
            }
        }

        CELibSNIIYieldsIntegrated[i].Metallicity = CELibSNIIYields(i,0).Metallicity;
        pCELibSNIIYieldsCurrentMetallicity = i;

        // Calculate Mr, Msw, E 

        CELibSNIIYieldsIntegrated[i].Mr = pCELibGetSNIIYieldsMrInGivenMassRange(MassL,MassU);
        CELibSNIIYieldsIntegrated[i].Msw = pCELibGetSNIIYieldsMswInGivenMassRange(MassL,MassU);
        CELibSNIIYieldsIntegrated[i].Msn = pCELibGetSNIIYieldsMsnInGivenMassRange(MassL,MassU);
        CELibSNIIYieldsIntegrated[i].Energy = pCELibGetSNIIYieldsEnergyInGivenMassRange(MassL,MassU);

        if(CELibRunParameters.TestMode){
            fprintf(stderr,"//\t Integrated SNII yields.\n");
            fprintf(stderr,"//\t\t Metallicity = %g\n",CELibSNIIYieldsIntegrated[i].Metallicity);
            fprintf(stderr,"//\t\t Mr = %g\n",CELibSNIIYieldsIntegrated[i].Mr);
            fprintf(stderr,"//\t\t Msw = %g\n",CELibSNIIYieldsIntegrated[i].Msw);
            fprintf(stderr,"//\t\t Msn = %g\n",CELibSNIIYieldsIntegrated[i].Msn);
            fprintf(stderr,"//\t\t Energy = %g\n",CELibSNIIYieldsIntegrated[i].Energy);
            fprintf(stderr,"//\t\t Mejecta = %g\n",CELibSNIIYieldsIntegrated[i].Msw+CELibSNIIYieldsIntegrated[i].Msn);
        }

        if((CELibRunParameters.SNIIYieldsTableID == CELibSNIIYieldsTableID_P98)&&(CELibRunParameters.SNIIYieldsModificationP98)){
            for(int j=0;j<CELibYield_Number;j++){
                pCELibSNIIYieldsCurrentElement = j;
                // CELibSNIIYieldsIntegrated[i].Elements[j] = 
                    // pCELibGetSNIIYieldsYieldMassInGivenMassRange(CELibRunParameters.SNIILowerMass,CELibRunParameters.SNIIUpperMass);
                CELibSNIIYieldsIntegrated[i].Elements[j] = 
                    pCELibGetSNIIYieldsYieldMassInGivenMassRange(MassL,MassU);
                if(j==CELibYield_C)  CELibSNIIYieldsIntegrated[i].Elements[j] *= 0.5;
                if(j==CELibYield_Mg) CELibSNIIYieldsIntegrated[i].Elements[j] *= 2.0;
                if(j==CELibYield_Fe) CELibSNIIYieldsIntegrated[i].Elements[j] *= 0.5;
            }
        } else {
            for(int j=0;j<CELibYield_Number;j++){
                pCELibSNIIYieldsCurrentElement = j;
                // CELibSNIIYieldsIntegrated[i].Elements[j] = 
                    // pCELibGetSNIIYieldsYieldMassInGivenMassRange(CELibRunParameters.SNIILowerMass,CELibRunParameters.SNIIUpperMass);
                CELibSNIIYieldsIntegrated[i].Elements[j] = 
                   pCELibGetSNIIYieldsYieldMassInGivenMassRange(MassL,MassU);
            }
        }

        if(CELibRunParameters.TestMode){
            for(int j=0;j<CELibYield_Number;j++){
                char Name[MaxCharactersInLine];
                CELibGetSNIIYieldElementName(j,Name);
                fprintf(stderr,"//\t\t %-2s %g\n",Name,CELibSNIIYieldsIntegrated[i].Elements[j]);
            }
            fprintf(stderr,"//\n");
        }
    }

    return ;
}


/*!
 * This function returns the integrated yields data in a given mass range.
 * Assume 1Msun SSP particle.
 */
struct CELibStructFeedbackOutput CELibGetSNIIYieldsIntegratedInGivenMassRange(const int MetalID, const double LowerMass_in_SolarMass, const double UpperMass_in_SolarMass){

    struct CELibStructFeedbackOutput SNIIFeedback;

    pCELibSNIIYieldsCurrentIMFType = CELibRunParameters.IMFType;

    double InitElements[CELibYield_Number];

    if(CELibRunParameters.SNIIYieldsTableID == CELibSNIIYieldsTableID_N13){
        if(CELibRunParameters.PopIIISNe == 1){
            if(CELibSNIIYields(MetalID,0).Metallicity <= CELibRunParameters.PopIIIMetallicity){
                if(CELibRunParameters.PopIIIIMF == 1){
                    pCELibSNIIYieldsCurrentIMFType = CELibIMF_Susa;
                    // MassU = CELibIMF[CELibIMF_Susa].MassMax;
                    CELibSetPrimordialMetallicity(1.0,InitElements); 
                } else {
                    pCELibSNIIYieldsCurrentIMFType = CELibRunParameters.IMFType;
                    // MassU = CELibIMF[CELibRunParameters.IMFType].MassMax;
                    CELibSetMetallicityWithSolarAbundancePattern(1.0,InitElements,CELibSNIIYields(MetalID,0).Metallicity); 
                }
            } else {
                pCELibSNIIYieldsCurrentIMFType = CELibRunParameters.IMFType;
                CELibSetMetallicityWithSolarAbundancePattern(1.0,InitElements,CELibSNIIYields(MetalID,0).Metallicity); 
            }
        }
    }

    pCELibSNIIYieldsCurrentMetallicity = MetalID;

    // double Mr = pCELibGetSNIIYieldsMrInGivenMassRange(LowerMass_in_SolarMass,UpperMass_in_SolarMass);
    double Msw = pCELibGetSNIIYieldsMswInGivenMassRange(LowerMass_in_SolarMass,UpperMass_in_SolarMass);
    double Msn = pCELibGetSNIIYieldsMsnInGivenMassRange(LowerMass_in_SolarMass,UpperMass_in_SolarMass);
    double Energy = pCELibGetSNIIYieldsEnergyInGivenMassRange(LowerMass_in_SolarMass,UpperMass_in_SolarMass);

    double Elements[CELibYield_Number];
    for(int j=0;j<CELibYield_Number;j++){
        Elements[j] = 0.e0;
    }
    if((CELibRunParameters.SNIIYieldsTableID == CELibSNIIYieldsTableID_P98)&&(CELibRunParameters.SNIIYieldsModificationP98)){
        for(int j=0;j<CELibYield_Number;j++){
            pCELibSNIIYieldsCurrentElement = j;
            Elements[j] = 
                pCELibGetSNIIYieldsYieldMassInGivenMassRange(LowerMass_in_SolarMass,
                        UpperMass_in_SolarMass);
            if(j==CELibYield_C)  Elements[j] *= 0.5;
            if(j==CELibYield_Mg) Elements[j] *= 2.0;
            if(j==CELibYield_Fe) Elements[j] *= 0.5;
        }
    } else {
        for(int j=0;j<CELibYield_Number;j++){
            pCELibSNIIYieldsCurrentElement = j;
            Elements[j] = 
                pCELibGetSNIIYieldsYieldMassInGivenMassRange(LowerMass_in_SolarMass,
                        UpperMass_in_SolarMass);
        }
    }


    // fprintf(stderr,"//\t\t Remnant mass = %g\n",Mr);
    // fprintf(stderr,"//\t\t SW mass = %g\n",Msw);
    // fprintf(stderr,"//\t\t SN mass = %g\n",Msn);
    // fprintf(stderr,"//\t\t Energy = %g\n",Energy);
    double EjectaMass = 0.e0;
    for(int i=0;i<CELibYield_Number;i++){
        char Name[MaxCharactersInLine];
        CELibGetSNIIYieldElementName(i,Name);
        Elements[i] += (Msw+Msn)*InitElements[i];
        EjectaMass += Elements[i];
        // fprintf(stderr,"//\t\t Elements %s = %g\n",Name,Elements[i]);
    }
    // fprintf(stderr,"//\t\t Ejecta mass = %g\n",EjectaMass);

    SNIIFeedback.Energy = Energy;
    SNIIFeedback.EjectaMass = EjectaMass;
    for(int i=0;i<CELibYield_Number;i++){
        SNIIFeedback.Elements[i] = Elements[i];
    }

    return SNIIFeedback;
}

struct CELibStructFeedbackOutput CELibGetSNIIYieldsIntegratedInGivenMassRangeZ(const double Metallicity, const double LowerMass_in_SolarMass, const double UpperMass_in_SolarMass){

    struct CELibStructFeedbackOutput SNIIFeedback;
    if(Metallicity <= CELibRunParameters.PopIIIMetallicity){
        SNIIFeedback = CELibGetSNIIYieldsIntegratedInGivenMassRange(0,LowerMass_in_SolarMass,UpperMass_in_SolarMass);
    } else if(Metallicity >= CELibSNIIYieldsZ[CELibSNIIYields_Metallicity-1]){
        SNIIFeedback = CELibGetSNIIYieldsIntegratedInGivenMassRange(CELibSNIIYields_Metallicity-1,
                            LowerMass_in_SolarMass,UpperMass_in_SolarMass);
    } else { // interpolation
        int IndexZ = CELibSNIIYields_Metallicity-2;
        for(int i=1;i<CELibSNIIYields_Metallicity;i++){
            if(CELibSNIIYieldsZ[i] > Metallicity){
                IndexZ = i-1;
                break;
            }
        }

        struct CELibStructFeedbackOutput SNIIFeedback1 = 
            CELibGetSNIIYieldsIntegratedInGivenMassRange(IndexZ,LowerMass_in_SolarMass,UpperMass_in_SolarMass);
        struct CELibStructFeedbackOutput SNIIFeedback2 = 
            CELibGetSNIIYieldsIntegratedInGivenMassRange(IndexZ+1,LowerMass_in_SolarMass,UpperMass_in_SolarMass);
        double EjectaMass = 0.e0;
        SNIIFeedback.Energy = (SNIIFeedback2.Energy-SNIIFeedback1.Energy)/(CELibSNIIYieldsZ[IndexZ+1]-CELibSNIIYieldsZ[IndexZ])
                                *(Metallicity-CELibSNIIYieldsZ[IndexZ])+SNIIFeedback1.Energy;
        for(int k=0;k<CELibYield_Number;k++){
            SNIIFeedback.Elements[k] = 
                (SNIIFeedback2.Elements[k]-SNIIFeedback1.Elements[k])/(CELibSNIIYieldsZ[IndexZ+1]-CELibSNIIYieldsZ[IndexZ])
                    *(Metallicity-CELibSNIIYieldsZ[IndexZ])+SNIIFeedback1.Elements[k];
            EjectaMass += SNIIFeedback.Elements[k];
        }
        SNIIFeedback.EjectaMass = EjectaMass;
        // gprint(SNIIFeedback1.EjectaMass);
        // gprint(SNIIFeedback2.EjectaMass);
        // gprint(SNIIFeedback.EjectaMass);
    }

    return SNIIFeedback;
}

/*!
 * This function writes the integrated yields data in files.
 */
static void pCElibWriteSNIIYieldsIntegratedTable(const char OutDir[]){

    pCELibSNIIYieldsCurrentIMFType = CELibRunParameters.IMFType;

    FILE *fp;
    char fname[MaxCharactersInLine];

    MakeDir(OutDir);

    for(int i=0;i<CELibSNIIYields_Metallicity;i++){
        if(CELibRunParameters.SNIIYieldsTableID == CELibSNIIYieldsTableID_P98){
            Snprintf(fname,"%s/SNIIYieldsIntegratedP98.%02d",OutDir,i);
        } else if(CELibRunParameters.SNIIYieldsTableID == CELibSNIIYieldsTableID_N13){
            Snprintf(fname,"%s/SNIIYieldsIntegratedN13.%02d",OutDir,i);
        }
        FileOpen(fp,fname,"w");
        fprintf(fp,"#Metallicity = %g\n",CELibSNIIYieldsIntegrated[i].Metallicity);
        fprintf(fp,"#Mr = %g\n",CELibSNIIYieldsIntegrated[i].Mr);
        fprintf(fp,"#Msw = %g\n",CELibSNIIYieldsIntegrated[i].Msw);
        fprintf(fp,"#Msn = %g\n",CELibSNIIYieldsIntegrated[i].Msn);
        fprintf(fp,"#Energy = %g\n",CELibSNIIYieldsIntegrated[i].Energy);
        fprintf(fp,"#Mejecta = %g\n",CELibSNIIYieldsIntegrated[i].Msw+CELibSNIIYieldsIntegrated[i].Msn);

        for(int j=0;j<CELibYield_Number;j++){
            char Name[MaxCharactersInLine];
            CELibGetSNIIYieldElementName(j,Name);
            fprintf(fp,"%g %s\n",CELibSNIIYieldsIntegrated[i].Elements[j],Name);
        }
        fclose(fp);
    }

    return ;
}


/*!
 * Write return mass fractions.
 */
static void pCElibWriteSNIIYieldsReturnMassFraction(const char OutDir[]){

    FILE *fp;
    char fname[MaxCharactersInLine];
    Snprintf(fname,"%s/SNIIRetrunMassFraction%d.dat",OutDir,CELibRunParameters.SNIIYieldsTableID);
    FileOpen(fp,fname,"w");

    for(int k=0;k<CELibIMF_NTypes;k++){
        pCELibSNIIYieldsCurrentIMFType = k;
        fprintf(fp,"IMF Type: %d\n",k);
        for(int i=0;i<CELibSNIIYields_Metallicity;i++){

            double Z = CELibSNIIYields(i,0).Metallicity;
            pCELibSNIIYieldsCurrentMetallicity = i;
            double Erf = pCELibGetSNIIYieldsMswInGivenMassRange(CELibRunParameters.SNIILowerMass,CELibRunParameters.SNIIUpperMass)
                    +pCELibGetSNIIYieldsMsnInGivenMassRange(CELibRunParameters.SNIILowerMass,CELibRunParameters.SNIIUpperMass);
            fprintf(fp,"Metallicity = %g\n",Z);
            fprintf(fp,"Erf = %g\n",Erf);
        }
        fprintf(fp,"\n");

    }
    fclose(fp);

    pCELibSNIIYieldsCurrentIMFType = CELibRunParameters.IMFType;

    return ;
}

/*!
 * Write hyper nova yields.
 */
static void pCElibWriteSNIIYieldsHyperNovaYields(const char OutDir[]){

    FILE *fp;
    char fname[MaxCharactersInLine];
    Snprintf(fname,"%s/SNIIHyperNova.dat",OutDir);
    FileOpen(fp,fname,"w");

    for(int i=0;i<CELibSNIIYields_Metallicity_N13HN;i++){
        int Nbin = 20;
        double MassMin = 10.0;
        double MassMax = 150.0;
        double dMass = (MassMax-MassMin)/Nbin;

        for(int k=0;k<Nbin;k++){
            double Mass = k*dMass + MassMin;
            fprintf(fp,"Mass = %g, Metallicity = %g\n",Mass,CELibSNIIYieldsN13HN[i][0].Metallicity);
            double Elements[CELibYield_Number];
            double Mr,Msw,Msn,Energy;
            int check =  pCELibGetSNIIYieldsHyperNovaYields(i,Mass,Elements,&Mr,&Msw,&Msn,&Energy);
            if(check == 1){
                fprintf(fp,"%d\n %g\n %g\n %g\n %g\n",check,Mr,Msw,Msn,Energy);
                for(int j=0;j<CELibYield_Number;j++){
                    fprintf(fp,"%g %s\n",Elements[j],SNIIYieldsElementName[j]);
                }
            }
        }
        fprintf(fp,"\n");

    }
    fclose(fp);

    return ;
}

static double pCELibSNIIYieldsMetallicityLower;
static double pCELibSNIIYieldsMetallicityUpper;
static double *pCELibSNIIYieldsMetallicityArray;
static double *pCELibSNIIYieldsNumberPerMassArray;
static double *pCELibSNIIYieldsMassStellarWindArray;
static double *pCELibSNIIYieldsMassSupernovaArray;

/*!
 * In this function, CELib copies frequently used values defined in the yields
 * structure to distinct arrays.
 */
static void pCELibCopySNIIYieldsRanges(void){

    pCELibSNIIYieldsMetallicityArray = realloc(pCELibSNIIYieldsMetallicityArray,sizeof(double)*CELibSNIIYields_Metallicity);
    pCELibSNIIYieldsNumberPerMassArray = realloc(pCELibSNIIYieldsNumberPerMassArray,sizeof(double)*CELibSNIIYields_Metallicity);
    pCELibSNIIYieldsMassStellarWindArray = realloc(pCELibSNIIYieldsMassStellarWindArray,sizeof(double)*CELibSNIIYields_Metallicity);
    pCELibSNIIYieldsMassSupernovaArray = realloc(pCELibSNIIYieldsMassSupernovaArray,sizeof(double)*CELibSNIIYields_Metallicity);

    if(CELibRunParameters.TestMode == 1){
        fprintf(stderr,"//\t Show CElib SNII Integrated info (fHN=%g)\n",CELibRunParameters.SNIIHyperNovaFraction);
    }
    for(int i=0;i<CELibSNIIYields_Metallicity;i++){
        pCELibSNIIYieldsMetallicityArray[i]     = CELibSNIIYieldsIntegrated[i].Metallicity;
        pCELibSNIIYieldsNumberPerMassArray[i]   = CELibSNIIYieldsIntegrated[i].Energy;
        pCELibSNIIYieldsMassStellarWindArray[i] = CELibSNIIYieldsIntegrated[i].Msw;
        pCELibSNIIYieldsMassSupernovaArray[i]   = CELibSNIIYieldsIntegrated[i].Msn;
        if(CELibRunParameters.TestMode == 1){
            fprintf(stderr,"//\t\t %d Z=%g, Energy=%g, SW = %g, SN = %g\n",
                    i,pCELibSNIIYieldsMetallicityArray[i],
                    pCELibSNIIYieldsNumberPerMassArray[i],
                    pCELibSNIIYieldsMassStellarWindArray[i],
                    pCELibSNIIYieldsMassSupernovaArray[i]);
        }
    }
    return ;
}


/*!
 * 
 */
static double pCELibGetSNIIYieldsStellarWindMass(double Metallicity){
    
    if(Metallicity <= pCELibSNIIYieldsMetallicityLower){
        return pCELibSNIIYieldsMassStellarWindArray[0];
    } else if(Metallicity >= pCELibSNIIYieldsMetallicityUpper){
        return pCELibSNIIYieldsMassStellarWindArray[CELibSNIIYields_Metallicity-1];
    } else {
        for(int i=1;i<CELibSNIIYields_Metallicity;i++){
            if(Metallicity > pCELibSNIIYieldsMetallicityArray[i]){
                return (pCELibSNIIYieldsMassStellarWindArray[i]-pCELibSNIIYieldsMassStellarWindArray[i-1])/
                    (pCELibSNIIYieldsMetallicityArray[i]-pCELibSNIIYieldsMetallicityArray[i-1])
                    *(Metallicity-pCELibSNIIYieldsMetallicityArray[i-1])+pCELibSNIIYieldsMassStellarWindArray[i-1];
            }
        }
    }
    return NONE;
}

/*!
 *
 */
static double pCELibGetSNIIYieldsSupernovaMass(double Metallicity){
    
    if(Metallicity <= pCELibSNIIYieldsMetallicityLower){
        return pCELibSNIIYieldsMassSupernovaArray[0];
    } else if(Metallicity >= pCELibSNIIYieldsMetallicityUpper){
        return pCELibSNIIYieldsMassSupernovaArray[CELibSNIIYields_Metallicity-1];
    } else {
        for(int i=1;i<CELibSNIIYields_Metallicity;i++){
            if(Metallicity > pCELibSNIIYieldsMetallicityArray[i]){
                return (pCELibSNIIYieldsMassSupernovaArray[i]-pCELibSNIIYieldsMassSupernovaArray[i-1])/
                    (pCELibSNIIYieldsMetallicityArray[i]-pCELibSNIIYieldsMetallicityArray[i-1])
                    *(Metallicity-pCELibSNIIYieldsMetallicityArray[i-1])+pCELibSNIIYieldsMassSupernovaArray[i-1];
            }
        }
    }
    return NONE;
}

/*!
 * This function returns energy, masses of individual elements, and the sum of
 * them (ejecta mass) released from SNeII.  It also returns the remnant mass.
 */
struct CELibStructFeedbackOutput CELibGetSNIIFeedback(struct CELibStructFeedbackInput Input){

    struct CELibStructFeedbackOutput SNIIFeedback;

    if(Input.Metallicity <= CELibSNIIYieldsIntegrated[0].Metallicity){
        int TableID = 0;
        if(Input.noPopIII){
            TableID ++;
        }
        SNIIFeedback.Energy = CELibSNIIYieldsIntegrated[TableID].Energy;
        for(int i=0;i<CELibYield_Number;i++)
            SNIIFeedback.Elements[i] = CELibSNIIYieldsIntegrated[TableID].Elements[i];
    } else if(Input.Metallicity >= CELibSNIIYieldsIntegrated[CELibSNIIYields_Metallicity-1].Metallicity){
        SNIIFeedback.Energy = CELibSNIIYieldsIntegrated[CELibSNIIYields_Metallicity-1].Energy;
        for(int i=0;i<CELibYield_Number;i++)
            SNIIFeedback.Elements[i] = CELibSNIIYieldsIntegrated[CELibSNIIYields_Metallicity-1].Elements[i];
    } else {
        int IndexMetallicity = 0;
        for(int i=1;i<CELibSNIIYields_Metallicity;i++){
            if(CELibSNIIYieldsIntegrated[i].Metallicity > Input.Metallicity){
                IndexMetallicity = i-1;
                break;
            }
        }
        double MetallicityEdges[] = {CELibSNIIYieldsIntegrated[IndexMetallicity].Metallicity, 
                                     CELibSNIIYieldsIntegrated[IndexMetallicity+1].Metallicity};
        for(int i=0;i<CELibYield_Number;i++){
            double GradYield = (CELibSNIIYieldsIntegrated[IndexMetallicity+1].Elements[i]-
                                        CELibSNIIYieldsIntegrated[IndexMetallicity].Elements[i])/
                                     (MetallicityEdges[1]-MetallicityEdges[0]);
            SNIIFeedback.Elements[i] = GradYield*(Input.Metallicity-MetallicityEdges[0]) + 
                                        CELibSNIIYieldsIntegrated[IndexMetallicity].Elements[i];
        }

        double GradEnergy = (CELibSNIIYieldsIntegrated[IndexMetallicity+1].Energy-
                                    CELibSNIIYieldsIntegrated[IndexMetallicity].Energy)/
                                 (MetallicityEdges[1]-MetallicityEdges[0]);
        SNIIFeedback.Energy = GradEnergy*(Input.Metallicity-MetallicityEdges[0]) + 
                                    CELibSNIIYieldsIntegrated[IndexMetallicity].Energy;
    }

    double Msw = pCELibGetSNIIYieldsStellarWindMass(Input.Metallicity);
    double Msn = pCELibGetSNIIYieldsSupernovaMass(Input.Metallicity);
    double Mstar = Input.Mass*Input.MassConversionFactor;
    SNIIFeedback.EjectaMass = 0.e0;
    for(int i=0;i<CELibYield_Number;i++){
        SNIIFeedback.Elements[i] *= Mstar;
        SNIIFeedback.Elements[i] += (Msw+Msn)*Mstar*(Input.Elements[i]/Input.Mass);
        SNIIFeedback.EjectaMass += SNIIFeedback.Elements[i];
    }
    SNIIFeedback.RemnantMass = Input.Mass*Input.MassConversionFactor - SNIIFeedback.EjectaMass;

    SNIIFeedback.Energy *= Input.Mass*Input.MassConversionFactor;

    return SNIIFeedback;
}

/*!
 * This function returns the total energy in units of erg released from a single
 * SSP particle whose mass is 1 Msun.
 */
double CELibGetSNIIIntegratedEnergy(const double Metallicity){

    if(Metallicity <= CELibSNIIYieldsIntegrated[0].Metallicity){
        return CELibSNIIYieldsIntegrated[0].Energy;
    } else if(Metallicity >= CELibSNIIYieldsIntegrated[CELibSNIIYields_Metallicity-1].Metallicity){
        return CELibSNIIYieldsIntegrated[CELibSNIIYields_Metallicity-1].Energy;
    } else {
        int IndexMetallicity = 0;
        for(int i=1;i<CELibSNIIYields_Metallicity;i++){
            if(CELibSNIIYieldsIntegrated[i].Metallicity > Metallicity){
                IndexMetallicity = i-1;
                break;
            }
        }
        double MetallicityEdges[] = {CELibSNIIYieldsIntegrated[IndexMetallicity].Metallicity, 
                                     CELibSNIIYieldsIntegrated[IndexMetallicity+1].Metallicity};
        double GradEnergy = (CELibSNIIYieldsIntegrated[IndexMetallicity+1].Energy-
                                    CELibSNIIYieldsIntegrated[IndexMetallicity].Energy)/
                                 (MetallicityEdges[1]-MetallicityEdges[0]);
        return GradEnergy*(Metallicity-MetallicityEdges[0]) + 
                            CELibSNIIYieldsIntegrated[IndexMetallicity].Energy;
    }
}


/*!
 * This function returns energy, masses of individual elements, and the sum of
 * them (ejecta mass) released from a single SNII.  It also returns the remnant
 * mass.
 */
struct CELibStructFeedbackStarbyStarOutput CELibGetSNIIFeedbackStarbyStar(struct CELibStructFeedbackStarbyStarInput Input){

    struct CELibStructFeedbackStarbyStarOutput SNIIStarbyStarFeedback;

    double Mstar = Input.Mass*Input.MassConversionFactor;

    if(Input.Metallicity <= CELibSNIIYieldsIntegrated[0].Metallicity){
        int TableID = 0;
        if(Input.noPopIII == 1){
            TableID ++;
        }

        double Msw = pCELibGetSNIIYieldsInterpolatedStellarWindMass(TableID,Mstar);
        double Msn = pCELibGetSNIIYieldsInterpolatedSuperNovaMass(TableID,Mstar);

        SNIIStarbyStarFeedback.EjectaMass = 0.e0;
        for(int i=0;i<CELibYield_Number;i++){
            SNIIStarbyStarFeedback.Elements[i] = pCELibGetSNIIYieldsInterpolatedYieldMass(TableID,i,Mstar)
                                      +(Msw+Msn)*(Input.Elements[i]/Input.Mass);
            SNIIStarbyStarFeedback.EjectaMass += SNIIStarbyStarFeedback.Elements[i];
        }
        SNIIStarbyStarFeedback.Energy = pCELibGetSNIIYieldsInterpolatedEnergy(TableID,Mstar);

    } else if(Input.Metallicity >= CELibSNIIYieldsIntegrated[CELibSNIIYields_Metallicity-1].Metallicity){
        double Msw = pCELibGetSNIIYieldsInterpolatedStellarWindMass(CELibSNIIYields_Metallicity-1,Mstar);
        double Msn = pCELibGetSNIIYieldsInterpolatedSuperNovaMass(CELibSNIIYields_Metallicity-1,Mstar);

        SNIIStarbyStarFeedback.EjectaMass = 0.e0;
        for(int i=0;i<CELibYield_Number;i++){
            SNIIStarbyStarFeedback.Elements[i] = pCELibGetSNIIYieldsInterpolatedYieldMass(CELibSNIIYields_Metallicity-1,i,Mstar)
                                      +(Msw+Msn)*(Input.Elements[i]/Input.Mass);
            SNIIStarbyStarFeedback.EjectaMass += SNIIStarbyStarFeedback.Elements[i];
        }
        SNIIStarbyStarFeedback.Energy = pCELibGetSNIIYieldsInterpolatedEnergy(CELibSNIIYields_Metallicity-1,Mstar);
    } else {
        int IndexMetallicity = 0;
        for(int i=1;i<CELibSNIIYields_Metallicity;i++){
            if(CELibSNIIYieldsIntegrated[i].Metallicity > Input.Metallicity){
                IndexMetallicity = i-1;
                break;
            }
        }

        double MetallicityEdges[] = {CELibSNIIYieldsIntegrated[IndexMetallicity].Metallicity, 
                                     CELibSNIIYieldsIntegrated[IndexMetallicity+1].Metallicity};

        double MswEdges[] = {pCELibGetSNIIYieldsInterpolatedStellarWindMass(IndexMetallicity,Mstar),
                             pCELibGetSNIIYieldsInterpolatedStellarWindMass(IndexMetallicity+1,Mstar)};
        double MsnEdges[] = {pCELibGetSNIIYieldsInterpolatedSuperNovaMass(IndexMetallicity,Mstar),
                             pCELibGetSNIIYieldsInterpolatedSuperNovaMass(IndexMetallicity+1,Mstar)};

        double GradMsw = (MswEdges[1]-MswEdges[0])/(MetallicityEdges[1]-MetallicityEdges[0]);
        double GradMsn = (MsnEdges[1]-MsnEdges[0])/(MetallicityEdges[1]-MetallicityEdges[0]);

        double Msw = GradMsw*(Input.Metallicity-MetallicityEdges[0])+MswEdges[0];
        double Msn = GradMsn*(Input.Metallicity-MetallicityEdges[0])+MsnEdges[0];


        SNIIStarbyStarFeedback.EjectaMass = 0.e0;
        for(int i=0;i<CELibYield_Number;i++){
            double GradYield = (pCELibGetSNIIYieldsInterpolatedYieldMass(IndexMetallicity+1,i,Mstar)-
                                        pCELibGetSNIIYieldsInterpolatedYieldMass(IndexMetallicity,i,Mstar))/
                                     (MetallicityEdges[1]-MetallicityEdges[0]);

            SNIIStarbyStarFeedback.Elements[i] = GradYield*(Input.Metallicity-MetallicityEdges[0]) + 
                                        pCELibGetSNIIYieldsInterpolatedYieldMass(IndexMetallicity,i,Mstar)
                                        + (Msw+Msn)*(Input.Elements[i]/Input.Mass);
            SNIIStarbyStarFeedback.EjectaMass += SNIIStarbyStarFeedback.Elements[i];
        }

        double GradEnergy = (pCELibGetSNIIYieldsInterpolatedEnergy(IndexMetallicity+1,Mstar)-
                                pCELibGetSNIIYieldsInterpolatedEnergy(IndexMetallicity,Mstar))/
                                 (MetallicityEdges[1]-MetallicityEdges[0]);
        SNIIStarbyStarFeedback.Energy = GradEnergy*(Input.Metallicity-MetallicityEdges[0]) + 
                                    pCELibGetSNIIYieldsInterpolatedEnergy(IndexMetallicity,Mstar);
    }

    SNIIStarbyStarFeedback.RemnantMass = Input.Mass*Input.MassConversionFactor - SNIIStarbyStarFeedback.EjectaMass;
    SNIIStarbyStarFeedback.Energy *= Input.Mass*Input.MassConversionFactor;

    return SNIIStarbyStarFeedback;
}
