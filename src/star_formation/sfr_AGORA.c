/*!
 * \copyright   This file is part of the public version of the AREPO code.
 * \copyright   Copyright (C) 2009-2019, Max-Planck Institute for Astrophysics
 * \copyright   Developed by Volker Springel (vspringel@MPA-Garching.MPG.DE) and
 *              contributing authors.
 * \copyright   Arepo is free software: you can redistribute it and/or modify
 *              it under the terms of the GNU General Public License as published by
 *              the Free Software Foundation, either version 3 of the License, or
 *              (at your option) any later version.
 *
 *              Arepo is distributed in the hope that it will be useful,
 *              but WITHOUT ANY WARRANTY; without even the implied warranty of
 *              MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *              GNU General Public License for more details.
 *
 *              A copy of the GNU General Public License is available under
 *              LICENSE as part of this program.  See also
 *              <https://www.gnu.org/licenses/>.
 *
 * \file        src/star_formation/sfr_AGORA.c
 * \date        05/2018
 * \brief       Star formation rate routines for the effective multi-phase
 *              model.
 * \details     contains functions:
 *                void cooling_and_starformation(void)
 *                double get_starformation_rate(int i)
 *                void set_units_sfr(void)
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#include "../gravity/forcetree.h"

#if defined(USE_SFR) && defined(AGORA_SF)

/*! \brief Main driver for star formation and gas cooling.
 *
 *  This function loops over all the active gas cells. If a given cell
 *  meets the criteria for star formation to be active the multi-phase
 *  model is activated, the properties of the cell are updated according to
 *  the latter and the star formation rate computed. In the other case, the
 *  standard isochoric cooling is applied to the gas cell by calling the
 *  function cool_cell() and the star formation rate is set to 0.
 *
 *  \return void
 */
void cooling_and_starformation(void)
{
  TIMER_START(CPU_COOLINGSFR);

  int idx, i, bin, flag;
  double dt, dtime, ne = 1;
  double unew, du;

  double dens;
    
  double t_freefall;  // Freefall timescale
  double sf_dens_threshold;  // Number density of neutral atomic hydrogen, code units - converted from parameter file value

  sf_dens_threshold=All.StarFormationNumberDensityThreshold*PROTONMASS/All.UnitDensity_in_cgs;
    
  /* note: assuming FULL ionization */
  double u_to_temp_fac =
      (4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC))) * PROTONMASS / BOLTZMANN * GAMMA_MINUS1 * All.UnitEnergy_in_cgs / All.UnitMass_in_g;

  /* clear the SFR stored in the active timebins */
  for(bin = 0; bin < TIMEBINS; bin++)
    if(TimeBinSynchronized[bin])
      TimeBinSfr[bin] = 0;

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].Mass == 0 && P[i].ID == 0)
        continue; /* skip cells that have been swallowed or eliminated */

      dens = SphP[i].Density;

      dt    = (P[i].TimeBinHydro ? (((integertime)1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;
      dtime = All.cf_atime * dt / All.cf_time_hubble_a;

      /* apply the temperature floor */

      unew = dmax(All.MinEgySpec, SphP[i].Utherm);

      if(unew < 0)
        terminate("Invalid Temperature: Task=%d i=%d unew=%g\n", ThisTask, i, unew);

      du = unew - SphP[i].Utherm;
      SphP[i].Utherm += du;
      SphP[i].Energy += All.cf_atime * All.cf_atime * du * P[i].Mass;

      cool_cell(i);


      /* check whether conditions for star formation are fulfilled.
       * f=1  normal cooling
       * f=0  star formation
       */

      flag = 1; /* default is normal cooling */

      /* enable star formation if gas is above SF density threshold */
      if(dens * All.cf_a3inv >= sf_dens_threshold)
        flag = 0;

      if(All.ComovingIntegrationOn)
        if(dens < All.OverDensThresh)
          flag = 1;

      if(P[i].Mass == 0) /* tracer particles don't form stars */
        flag = 1;

      if(flag == 1)
        SphP[i].Sfr = 0;

      /* active star formation */
      if(flag == 0)
        {
          SphP[i].Ne = (HYDROGEN_MASSFRAC + 1) / 2 / HYDROGEN_MASSFRAC; /* note: assuming FULL ionization */

          if(dt > 0)
            {
              if(P[i].TimeBinHydro) /* upon start-up, we need to protect against dt==0 */
                {
                  unew = SphP[i].Utherm;

                  du = unew - SphP[i].Utherm;
                  if(unew < All.MinEgySpec)
                    du = All.MinEgySpec - SphP[i].Utherm;

                  SphP[i].Utherm += du;
                  SphP[i].Energy += All.cf_atime * All.cf_atime * du * P[i].Mass;

#ifdef OUTPUT_COOLHEAT
                  if(dtime > 0)
                    SphP[i].CoolHeat = du * P[i].Mass / dtime;
#endif /* #ifdef OUTPUT_COOLHEAT */

                  set_pressure_of_cell(i);
                }
            }

            //  Compute the freefall time
            t_freefall=sqrt(3.*M_PI/32/All.G/SphP[i].Density); // freefall time in code units

            SphP[i].Sfr=All.StarFormationEfficiency*P[i].Mass/t_freefall * (All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);
            
            TimeBinSfr[P[i].TimeBinHydro] += SphP[i].Sfr;
        }
    } /* end of main loop over active particles */

  TIMER_STOP(CPU_COOLINGSFR);
}

/*! \brief Return the star formation rate associated with the gas cell i.
 *
 *  \param[in] i the index of the gas cell.
 *
 *  \return star formation rate in solar masses / yr.
 */
double get_starformation_rate(int i)
{
  if(RestartFlag == 3)
    return SphP[i].Sfr;

  double rateOfSF;
  int flag;
  double tsfr;
  /* note: assuming FULL ionization */
  double u_to_temp_fac =
      (4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC))) * PROTONMASS / BOLTZMANN * GAMMA_MINUS1 * All.UnitEnergy_in_cgs / All.UnitMass_in_g;

  double t_freefall;  // Freefall timescale
  double sf_dens_threshold;  // Number density of neutral atomic hydrogen, code units - converted from parameter file value

  sf_dens_threshold=All.StarFormationNumberDensityThreshold*PROTONMASS*pow(All.UnitDensity_in_cgs,-1.);

  flag   = 1; /* default is normal cooling */
    
  if(SphP[i].Density * All.cf_a3inv >= sf_dens_threshold)
    flag = 0;

  if(All.ComovingIntegrationOn)
    if(SphP[i].Density < All.OverDensThresh)
      flag = 1;

  if(flag == 1)
    return 0;

  t_freefall=sqrt(3.*M_PI/32/All.G/SphP[i].Density); // freefall time in code units

  rateOfSF=All.StarFormationEfficiency*P[i].Mass/t_freefall;

  /* convert to solar masses per yr */
  rateOfSF *= (All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);

  return rateOfSF;
}

/*! \brief Set the appropriate units for the parameters of the multi-phase
 *         model.
 *
 *  \return void
 */
void set_units_sfr(void)
{
  double meanweight;

  All.OverDensThresh = All.CritOverDensity * All.OmegaBaryon * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);

  All.PhysDensThresh = All.CritPhysDensity * PROTONMASS / HYDROGEN_MASSFRAC / All.UnitDensity_in_cgs;
}

#endif /* #if defined(USE_SFR) && defined(AGORA_SF) */
