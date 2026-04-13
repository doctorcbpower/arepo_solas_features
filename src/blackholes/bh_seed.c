#include <gsl/gsl_math.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#include "../domain/domain.h"

#include "../fof/fof.h"

/*! \brief Spawn a black hole particle from the densest gas cell in a FOF group.
 *
 *  Finds the densest gas cell belonging to this group across all MPI tasks,
 *  then the owning task spawns a Type 5 BH particle from it, conserving mass.
 *
 *  \param[in] grp_index  Index into Group[] for this group.
 *  \return void
 */

#ifdef BLACKHOLE_SEEDING
void seed_black_hole_in_group(int grp_index, int *n_seeded)
{
  /* grp_index == -1 means this task has no group to seed but must
   * still participate in the MPI_Allreduce collective */
  
  MyIDType target_minid = (grp_index >= 0) ? Group[grp_index].MinID : 0;

  double local_max_density = -1.0;
  int    local_best_index  = -1;

  if(grp_index >= 0)
    {
      for(int i = 0; i < NumPart; i++)
        {
          if(P[i].Type != 0) continue;
          if(P[i].Mass == 0 && P[i].ID == 0) continue;
          if(FOF_PList[i].MinID != target_minid) continue;
          if(SphP[i].Density > local_max_density)
            {
              local_max_density = SphP[i].Density;
              local_best_index  = i;
            }
        }
    }

  struct { double density; int task; } local_val, global_val;
  local_val.density = (local_best_index >= 0) ? local_max_density : -1.0;
  local_val.task    = ThisTask;

  MPI_Allreduce(&local_val, &global_val, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);

  /* Rest of spawning logic unchanged - only winning task acts */
  if(grp_index < 0 || global_val.density < 0 || ThisTask != global_val.task)
    {
      mpi_printf("FOF_SEEDING: WARNING - no gas cell found for group MinID=%llu, skipping seed.\n",
                 (unsigned long long)target_minid);
      return;
    }

  /* --- Phase 3: winning task spawns the BH --- */
  if(ThisTask != global_val.task)
    return;  /* nothing to do on non-winning tasks */

  int igas = local_best_index;

  /* Check we have room — same pattern as make_star() */
  if(NumPart >= All.MaxPart)
    terminate("FOF_SEEDING: no space to spawn BH (NumPart=%d MaxPart=%d)", NumPart, All.MaxPart);

#ifdef BLACKHOLES
  /* Check BhP array has room */
  if(NumBhs >= All.MaxPartBhs)
    {
      All.MaxPartBhs = (int)(1.25 * All.MaxPartBhs) + 1;
      reallocate_memory_maxpartbhs();
    }
#endif

  int ibh = NumPart;   /* index of new BH in P[] */

  /* Copy base particle data from gas cell */
  P[ibh] = P[igas];
  P[ibh].Type          = 5;
  P[ibh].SofteningType = All.SofteningTypeOfPartType[5];
  P[ibh].Mass          = All.BlackHoleSeedMass;

#ifdef INDIVIDUAL_GRAVITY_SOFTENING
  if(((1 << P[ibh].Type) & (INDIVIDUAL_GRAVITY_SOFTENING)))
    P[ibh].SofteningType = get_softening_type_from_mass(P[ibh].Mass);
#endif

  /* Give it a unique ID */
  if(All.MaxID == 0)
    calculate_maxid();
  All.MaxID++;
  P[ibh].ID = All.MaxID;

  /* Register in gravity timebin — same as spawn_star_from_cell() */
  timebin_add_particle(&TimeBinsGravity, ibh, igas, P[ibh].TimeBinGrav,
                       TimeBinSynchronized[P[ibh].TimeBinGrav]);

  /* Reduce gas cell mass conservatively */
  double fac = (P[igas].Mass - All.BlackHoleSeedMass) / P[igas].Mass;
  P[igas].Mass          *= fac;
  SphP[igas].Energy     *= fac;
  SphP[igas].Momentum[0] *= fac;
  SphP[igas].Momentum[1] *= fac;
  SphP[igas].Momentum[2] *= fac;
#ifdef METALS
  SphP[igas].Metals *= fac;
#endif
#ifdef MAXSCALARS
  for(int s = 0; s < N_Scalar; s++)
    *(MyFloat *)(((char *)(&SphP[igas])) + scalar_elements[s].offset_mass) *= fac;
#endif

#ifdef BLACKHOLES
  /* Initialise BhP entry */
  P[ibh].BhID             = NumBhs;
  BhP[NumBhs].PID         = ibh;
  BhP[NumBhs].Hsml        = cbrt((3.0 * SphP[igas].Volume) / (4.0 * M_PI));
  BhP[NumBhs].Density      = SphP[igas].Density;
  BhP[NumBhs].NgbMass      = 0.0;
  BhP[NumBhs].NgbMassFeed  = 0.0;
  BhP[NumBhs].DensityFlag  = 0;
  BhP[NumBhs].TimeBinBh    = 0;
#ifdef BONDI_ACCRETION
  BhP[NumBhs].AccretionRate = 0.0;
  BhP[NumBhs].MassToDrain   = 0.0;
  for(int k = 0; k < 3; k++)
    {
      BhP[NumBhs].VelocityGas[k]         = 0.0;
      BhP[NumBhs].VelocityGasCircular[k] = 0.0;
      BhP[NumBhs].AngularMomentum[k]     = 0.0;
    }
  BhP[NumBhs].InternalEnergyGas = 0.0;
#endif
#ifdef INFALL_ACCRETION
  BhP[NumBhs].Accretion = 0.0;
#endif
  NumBhs++;
#endif /* BLACKHOLES */

  NumPart++;
  // All.TotNumPart++;
#ifdef BLACKHOLES
  // All.TotNumBhs++;
#endif

  mpi_printf("FOF_SEEDING: Seeded BH (ID=%llu) in group MinID=%llu, "
             "mass=%g at task %d\n",
             (unsigned long long)P[ibh].ID,
             (unsigned long long)target_minid,
             All.BlackHoleSeedMass, ThisTask);


  (*n_seeded)++;
}
#endif