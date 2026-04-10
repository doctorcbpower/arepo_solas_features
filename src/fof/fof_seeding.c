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
 * \file        src/fof/fof_seeding.c
 * \date        20/10/2025
 * \brief       Parallel friend of friends (FoF) group finder.
 * \details     contains functions:
 *                void fof_seeding(int num)
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 24.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include "fof.h"

#include <gsl/gsl_math.h>
#include <inttypes.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "../domain/domain.h"
#include "../main/allvars.h"
#include "../main/proto.h"
#include "../subfind/subfind.h"

#if defined(HALO_SEEDING) && defined(FOF)

#ifdef BLACKHOLE_SEEDING
  #ifndef BLACKHOLES
    #define BLACKHOLES // Temporary safeguard for the header
  #endif
  #include "../blackholes/bh_proto.h"
#endif

  int is_halo_seeded(MyIDType minid)
{
    // binary search
    int lo = 0, hi = All.NSeededHalos - 1;
    while(lo <= hi) {
        int mid = (lo + hi) / 2;
        if(SeededHaloIDs[mid] == minid) return 1;
        else if(SeededHaloIDs[mid] < minid) lo = mid + 1;
        else hi = mid - 1;
    }
    return 0;
}

void mark_halo_seeded(MyIDType minid)
{
    // Grow array if needed
    if(All.NSeededHalos == All.MaxSeededHalos) {
        All.MaxSeededHalos = 2 * All.MaxSeededHalos + 64;
        SeededHaloIDs = myrealloc(SeededHaloIDs, All.MaxSeededHalos * sizeof(MyIDType));
    }
    // Insert in sorted position
    int pos = All.NSeededHalos;
    while(pos > 0 && SeededHaloIDs[pos-1] > minid)  {
        SeededHaloIDs[pos] = SeededHaloIDs[pos-1];
        pos--;
    }
    SeededHaloIDs[pos] = minid;
    All.NSeededHalos++;
}

static MyIDType *MinID;
static int *Head, *Len, *Next, *Tail, *MinIDTask;

/*! \brief Main routine to execute the friend of friends group finder.
 *
 *  If called with num == -1 as argument, only FOF is carried out and no group
 *  catalogs are saved to disk. If num >= 0, the code will store the
 *  group/subgroup catalogs, and bring the particles into output order.
 *  In this case, the calling routine (which is normally savepositions()) will
 *  need to free PS[] and bring the particles back into the original order,
 *  as well as reestablished the mesh.
 *
 *  \param[in] num Index of output; if negative, no output written.
 *
 *  \return void
 */
void fof_seeding(void)
{
  int i, start, lenloc, largestgroup;
  double t0, t1, cputime;

  TIMER_START(CPU_FOF);

  mpi_printf("FOF_SEEDING: Begin to compute FoF group catalogue...  (presently allocated=%g MB)\n", AllocatedBytes / (1024.0 * 1024.0));

  ngb_treefree();

  domain_free();

  domain_Decomposition();

  ngb_treeallocate();
  ngb_treebuild(NumGas);

  /* check */
  for(i = 0; i < NumPart; i++)
    if((P[i].Mass == 0 && P[i].ID == 0) || (P[i].Type == 4 && P[i].Mass == 0) || (P[i].Type == 5 && P[i].Mass == 0))
      terminate("this should not happen");

  /* this structure will hold auxiliary information for each particle, needed only during group finding */
  PS = (struct subfind_data *)mymalloc_movable(&PS, "PS", All.MaxPart * sizeof(struct subfind_data));

  memset(PS, 0, NumPart * sizeof(struct subfind_data));

  /* First, we save the original location of the particles, in order to be able to revert to this layout later on */
  for(i = 0; i < NumPart; i++)
    {
      PS[i].OriginTask  = ThisTask;
      PS[i].OriginIndex = i;
    }

  fof_OldMaxPart    = All.MaxPart;
  fof_OldMaxPartSph = All.MaxPartSph;

  LinkL = fof_get_comoving_linking_length(); // in fof.c

  mpi_printf("FOF_SEEDING: Comoving linking length: %g    (presently allocated=%g MB)\n", LinkL, AllocatedBytes / (1024.0 * 1024.0));

  MinID     = (MyIDType *)mymalloc("MinID", NumPart * sizeof(MyIDType));
  MinIDTask = (int *)mymalloc("MinIDTask", NumPart * sizeof(int));

  Head = (int *)mymalloc("Head", NumPart * sizeof(int));
  Len  = (int *)mymalloc("Len", NumPart * sizeof(int));
  Next = (int *)mymalloc("Next", NumPart * sizeof(int));
  Tail = (int *)mymalloc("Tail", NumPart * sizeof(int));

#ifdef HIERARCHICAL_GRAVITY
  timebin_make_list_of_active_particles_up_to_timebin(&TimeBinsGravity, All.HighestOccupiedTimeBin);
#endif /* #ifdef HIERARCHICAL_GRAVITY */

  construct_forcetree(0, 0, 1, All.HighestOccupiedTimeBin); /* build tree for all particles */
//**
//#if defined(SUBFIND)
//  subfind_density_hsml_guess();
//#endif /* #if defined(SUBFIND) */

  /* initialize link-lists */
  for(i = 0; i < NumPart; i++)
    {
      Head[i] = Tail[i] = i;
      Len[i]            = 1;
      Next[i]           = -1;
      MinID[i]          = P[i].ID;
      MinIDTask[i]      = ThisTask;
    }

  /* call routine to find primary groups */
  cputime = fof_find_groups(MinID, Head, Len, Next, Tail, MinIDTask);
  mpi_printf("FOF_SEEDING: group finding took = %g sec\n", cputime);

//#ifdef FOF_SECONDARY_LINK_TARGET_TYPES
//  myfree(Father);
//  myfree(Nextnode);
//  myfree(Tree_Points);
//
//  /* now rebuild the tree with all the types selected as secondary link targets */
//  construct_forcetree(0, 0, 2, All.HighestOccupiedTimeBin);
//#endif /* #ifdef FOF_SECONDARY_LINK_TARGET_TYPES */

//#ifdef HIERARCHICAL_GRAVITY
//  timebin_make_list_of_active_particles_up_to_timebin(&TimeBinsGravity, All.HighestActiveTimeBin);
//#endif /* #ifdef HIERARCHICAL_GRAVITY */
//
//  /* call routine to attach secondary particles/cells to primary groups */
//  cputime = fof_find_nearest_dmparticle(MinID, Head, Len, Next, Tail, MinIDTask);
//
//  mpi_printf("FOF_SEEDING: attaching gas and star particles to nearest dm particles took = %g sec\n", cputime);

  myfree(Father);
  myfree(Nextnode);
  myfree(Tree_Points);

  force_treefree();

  myfree(Tail);
  myfree(Next);
  myfree(Len);
//  myfree(Head);
//  myfree(MinIDTask);
//  myfree(MinID);

  t0 = second();

  FOF_PList = (struct fof_particle_list *)mymalloc_movable(&FOF_PList, "FOF_PList", NumPart * sizeof(struct fof_particle_list));

  for(i = 0; i < NumPart; i++)
    {
      FOF_PList[i].MinID     = MinID[Head[i]];
      FOF_PList[i].MinIDTask = MinIDTask[Head[i]];
      FOF_PList[i].Pindex    = i;
    }

  myfree_movable(Head);
  myfree_movable(MinIDTask);
  myfree_movable(MinID);
    
  FOF_GList = (struct fof_group_list *)mymalloc_movable(&FOF_GList, "FOF_GList", sizeof(struct fof_group_list) * NumPart);

  fof_compile_catalogue(); // in fof.c

  t1 = second();
  mpi_printf("FOF_SEEDING: compiling local group data and catalogue took = %g sec\n", timediff(t0, t1));

  MPI_Allreduce(&Ngroups, &TotNgroups, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  sumup_large_ints(1, &Nids, &TotNids);

  if(TotNgroups > 0)
    {
      int largestloc = 0;

      for(i = 0; i < NgroupsExt; i++)
        if(FOF_GList[i].LocCount + FOF_GList[i].ExtCount > largestloc)
          largestloc = FOF_GList[i].LocCount + FOF_GList[i].ExtCount;
      MPI_Allreduce(&largestloc, &largestgroup, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    }
  else
    largestgroup = 0;

  mpi_printf("FOF_SEEDING: Total number of FOF groups with at least %d particles: %d\n", FOF_GROUP_MIN_LEN, TotNgroups);
  mpi_printf("FOF_SEEDING: Largest FOF group has %d particles.\n", largestgroup);
  mpi_printf("FOF_SEEDING: Total number of particles in FOF groups: %lld\n", TotNids);

  t0 = second();

  MaxNgroups = 2 * imax(NgroupsExt, TotNgroups / NTask + 1);

  Group = (struct group_properties *) mymalloc_movable(&Group, "Group", sizeof(struct group_properties) * MaxNgroups);

  mpi_printf("FOF_SEEDING: group properties are now allocated.. (presently allocated=%g MB)\n", AllocatedBytes / (1024.0 * 1024.0));

  for(i = 0, start = 0; i < NgroupsExt; i++)
    {
      while(FOF_PList[start].MinID < FOF_GList[i].MinID)
        {
          start++;
          if(start > NumPart)
            terminate("start > NumPart");
        }

      if(FOF_PList[start].MinID != FOF_GList[i].MinID)
        terminate("ID mismatch");

      for(lenloc = 0; start + lenloc < NumPart;)
        if(FOF_PList[start + lenloc].MinID == FOF_GList[i].MinID)
          lenloc++;
        else
          break;

      Group[i].MinID     = FOF_GList[i].MinID;
      Group[i].MinIDTask = FOF_GList[i].MinIDTask;

      fof_compute_group_properties(i, start, lenloc);

      start += lenloc;
    }

  fof_exchange_group_data();

  fof_finish_group_properties();

  t1 = second();
  mpi_printf("FOF_SEEDING: computation of group properties took = %g sec\n", timediff(t0, t1));
//
//  fof_assign_group_numbers();
//
  mpi_printf("FOF_SEEDING: Finished computing FoF groups.  (presently allocated=%g MB)\n", AllocatedBytes / (1024.0 * 1024.0));
//
#ifdef BLACKHOLE_SEEDING

  for(int n=0; n<Ngroups;n++)
      fprintf(stderr,"Group mass[%d]: %d (ThisTask:%d) \n",n,Group[n].Len,ThisTask);

    for(int n = 0; n < Ngroups; n++)
    {
      if(Group[n].Mass < All.MinHaloMassForBlackHoleSeeding)
        continue;
      if(is_halo_seeded(Group[n].MinID))
        continue;
      seed_black_hole_in_group(n);
      mark_halo_seeded(Group[n].MinID);

  }
#endif  

  myfree_movable(FOF_GList);
  myfree_movable(FOF_PList);

//#ifdef SUBFIND
//  TIMER_STOP(CPU_FOF);
//
//  subfind(0);
//
//  TIMER_START(CPU_FOF);
//#else  /* #ifdef SUBFIND */
//**  Nsubgroups    = 0;
//**  TotNsubgroups = 0;

  myfree_movable(Group);   
  myfree(PS);              

  TIMER_STOP(CPU_FOF);

  mpi_printf("FOF_SEEDING: All FOF related work finished...\n");
}
/*! \brief Calculate dynamical time at a given epoch.
 *
 *  \return Increment in expansion factor
 */
double fof_seeding_get_time_increment(void)
{
    double hubble=hubble_function(All.Time);
    return 2./hubble/sqrtf(200.);
}

#endif /* of FOF */
