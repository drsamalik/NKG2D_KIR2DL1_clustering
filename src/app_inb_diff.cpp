/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPPARKS directory.
------------------------------------------------------------------------- */
#include "math.h"
#include "mpi.h"
#include "stdlib.h"
#include "string.h"
#include "app_inb_diff.h"
#include "solve.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"
#include <ctime>
using namespace SPPARKS_NS;

enum{LOCAL,NBR};

#define DELTAEVENT 100000

/* ---------------------------------------------------------------------- */

AppInbDiff::AppInbDiff(SPPARKS *spk, int narg, char **arg) : 
  AppLattice(spk,narg,arg)
{
  ninteger = MAX_SPECIES;
  ndouble = 0;
  delpropensity = 1;
  delevent = 1;
  allow_kmc = 1;
  allow_rejection = 0;

  //total_pvav = 0;

  create_arrays();

  firsttime = 1;
  esites = NULL;
  echeck = NULL;
  events = NULL;
  ir_cluster_center = NULL;
  nevents = 0;
  maxevent = 0;
  firstevent = NULL;

  // species lists
  nspecies = 0;
  sname = NULL;

  // reaction lists
  nreactions = 0;
  rname = NULL;
  localReactants = NULL;
  nbrReactants = NULL;
  localDeltaPop = NULL;
  nbrDeltaPop = NULL;
  rate = NULL;
  rxnStyle = NULL;

  rxn_count = NULL; // count of how many times each reaction takes place

  // parse arguments for InbDiff class only, not children
  if (strcmp(style,"inb/diff") != 0) return;
  if (narg != 1) error->all(FLERR,"Illegal app_style command");
}

/* ---------------------------------------------------------------------- */

AppInbDiff::~AppInbDiff()
{
  delete [] esites;
  delete [] echeck;
  delete [] ir_cluster_center;
  memory->sfree(events);
  memory->destroy(firstevent);

  memory->destroy(sname);
  memory->destroy(rname);
  memory->destroy(localReactants);
  memory->destroy(nbrReactants);
  memory->destroy(localDeltaPop);
  memory->destroy(nbrDeltaPop);
  memory->destroy(rate);
  memory->destroy(rxnStyle);
  memory->destroy(rxn_count);
}

/* ---------------------------------------------------------------------- */

void AppInbDiff::input_app(char *command, int narg, char **arg)
{
  if (strcmp(command,"add_rxn") == 0) add_rxn(narg,arg);
  else if (strcmp(command,"add_species") == 0) add_species(narg,arg);
  else error->all(FLERR,"Unrecognized command");
}

/* ----------------------------------------------------------------------
   set site value ptrs each time iarray/darray are reallocated
------------------------------------------------------------------------- */

void AppInbDiff::grow_app()
{
  population = iarray;
  total_pvav = 0;
}

/* ----------------------------------------------------------------------
   initialize before each run
   check validity of site values
------------------------------------------------------------------------- */

void AppInbDiff::init_app()
{
  if (firsttime) {
    firsttime = 0;

    echeck = new int[nlocal];
    memory->create(firstevent,nlocal,"app:firstevent");

    // esites must be large enough for i, it's 1st neighbors and their 1st neighbors
    esites = new int[2 + 2*maxneigh]; //personalized
  }

  // zero reaction counts
  delete [] rxn_count;
  memory->create(rxn_count,nreactions,"inb/diff:rxn_count");
  for (int m = 0; m < nreactions; m++) rxn_count[m] = 0;

  // site validity
  int flag = 0;
  for (int ispecies = 0; ispecies < MAX_SPECIES; ispecies++) {
    for (int i = 0; i < nlocal; i++) {
      if (population[ispecies][i] < 0) flag = 1;
    }
  }
  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall) error->all(FLERR,"One or more sites have invalid values");
//  for (int nu = 0; nu <2312; nu++ ) {
//    total_pvav = total_pvav + population[11][nu]+population[14][nu]+population[24][nu]+population[25][nu]+ population[31][nu]+population[32][nu]+ population[33][nu]+population[34][nu];
//  }
/*  std::cout << "\n";
  std::cout << "sname_34=";
  std::cout << sname[34]; 
  std::cout << "\n"; 
  std::cout << "rname_34=";
  std::cout << rname[1]; 
  std::cout << "\n"; 
  std::cout << "rxn_style_34=";
  std::cout << rxnStyle[54]; 
  std::cout << "\n"; 
  std::cout << "delts_pop_l=";
  std::cout << localDeltaPop[54][0]; 
  std::cout << "\n"; 
  std::cout << "delts_pop_a=";
  std::cout << localDeltaPop[54][5]; 
  std::cout << "\n"; 
  std::cout << "delts_pop_r=";
  std::cout << localDeltaPop[54][28]; 
  std::cout << "\n"; 
  std::cout << "delts_pop_k=";
  std::cout << localDeltaPop[54][8]; 
  std::cout << "\n"; 
  std::cout << "delts_pop_D=";
  std::cout << localDeltaPop[54][30]; 
  std::cout << "\n"; */
}

/* ---------------------------------------------------------------------- */

void AppInbDiff::setup_app()
{
  for (int i = 0; i < nlocal; i++) echeck[i] = 0;

  // clear event list

  nevents = 0;
  for (int i = 0; i < nlocal; i++) firstevent[i] = -1;
  for (int i = 0; i < maxevent; i++) events[i].next = i+1;
  freeevent = 0;
  
  //int lattice_size_with_wall = 17*2;
  //int lattice_size = lattice_size_with_wall-4;
  //srand(0);
  
  /*int total_cluster_num = 20;
  int min_voxel_ir_cluster = (lattice_size_with_wall * 3) + 4; 
  int max_voxel_ir_cluster = (lattice_size_with_wall * 4) - 3; 
  int range_voxel = max_voxel_ir_cluster - min_voxel_ir_cluster;

  ir_cluster_center = new int[total_cluster_num];
  int *ir_choose = new int[range_voxel+1];
  int *ir_choose1 = new int[range_voxel+1];
  for (int kk=0; kk<(range_voxel+1); kk++) {
    ir_choose[kk] = min_voxel_ir_cluster + (kk * lattice_size_with_wall);
    ir_choose1[kk] = max_voxel_ir_cluster  + (kk * lattice_size_with_wall);
  }

//  ir_cluster_center[0] = 50; //45 + rand()%8;    
//  ir_cluster_center[1] = 150; // 145 + rand()%8;    
  srand(0);
  for (int kk=0; kk< (total_cluster_num); kk++) {
    if(kk < (total_cluster_num/4)) {
      ir_cluster_center[kk] = min_voxel_ir_cluster + rand()%(range_voxel + 1);    
      int ll =0;
    } else if (kk >= (total_cluster_num/4) && kk < (total_cluster_num/4)*2 ) {
      ir_cluster_center[kk] = (min_voxel_ir_cluster + (lattice_size_with_wall * range_voxel)) + rand()%(range_voxel + 1);    
      int ll =0;
    } else if (kk >= (total_cluster_num/4)*2 && kk < (total_cluster_num/4)*3 ) {    
      int r11 = rand() % (range_voxel+1); 
      ir_cluster_center[kk] = ir_choose[r11];
      int ll =0;
    } else if (kk >= (total_cluster_num/4)*3 && kk < (total_cluster_num/4)*4 ) {    
      int r22 = rand() % (range_voxel+1); 
      ir_cluster_center[kk] = ir_choose1[r22];         
      int ll =0;
    }
  }
  delete [] ir_choose;
  delete [] ir_choose1;*/

}

/* ----------------------------------------------------------------------
   compute energy of site
------------------------------------------------------------------------- */

double AppInbDiff::site_energy(int i)
{
  return 0.0;
}

/* ----------------------------------------------------------------------
   KMC method
   compute total propensity of owned site summed over possible events
------------------------------------------------------------------------- */

double AppInbDiff::site_propensity(int i)
{
  int j,k,m;

  // valid single, double, triple events are in tabulated lists
  // propensity for each event is input by user

  clear_events(i);

  double proball = 0.0;
  double my_rate, multiplicity;
  int ispecies, flag;

  // populate event list for this site
  for (m = 0; m < nreactions; m++) {

    // Local reactions
    if (rxnStyle[m] == LOCAL){
      multiplicity = 1.0; 
      flag = 0;

      // check that reactants are present
      for (ispecies=0; ispecies < MAX_SPECIES; ispecies++) {
        if (localReactants[m][ispecies]) { 
          if (population[ispecies][i] < localReactants[m][ispecies]){
            flag = 1;
            break;
          }

          // calculate multiplicity:
          // population[ispecies]i[i] choose localReactants[m][ispecies]
          for (k = 1; k<=localReactants[m][ispecies]; k++)
            multiplicity *= (double)(population[ispecies][i]+1-k)/(double)k;
        }
      }
      if (flag) continue;

      my_rate = rate[m] * multiplicity * custom_multiplier(i,LOCAL,m,-1);
      if (my_rate == 0.0) continue;
      add_event(i, LOCAL, m, my_rate, -1);
      proball += my_rate;

    // Neighbor reactions
    } else { // rxnStyle[m] == NBR
      for (int jj = 0; jj < numneigh[i]; jj++) {
        j = neighbor[i][jj];
        multiplicity = 1.0; 
        flag = 0;

        // check that reactants are present
        for (ispecies=0; ispecies < MAX_SPECIES; ispecies++) {
          // local
          if (localReactants[m][ispecies]) {
            if (population[ispecies][i] < localReactants[m][ispecies]){
              flag = 1;
              break;
            }

            // calculate multiplicity:
            // population[ispecies][i] choose localReactants[m][ispecies]
            for (k = 1; k<=localReactants[m][ispecies]; k++)
              multiplicity *= (double)(population[ispecies][i]+1-k)/(double)k;
          }

          // neighbor
          if (nbrReactants[m][ispecies]) {
            if (population[ispecies][j] < nbrReactants[m][ispecies]){
              flag = 1;
              break;
            }

            // calculate multiplicity:
            // population[ispecies][j] choose nbrReactants[m][ispecies]
            for (k = 1; k<=nbrReactants[m][ispecies]; k++)
              multiplicity *= (double)(population[ispecies][j]+1-k)/(double)k;
          }
        }
        if (flag) continue;
        my_rate = rate[m] * multiplicity * custom_multiplier(i,NBR,m,j);
        if (my_rate == 0.0) continue;
        add_event(i, NBR, m, my_rate, j);
        proball += my_rate;
      }
    }
  }

  return proball;
}

/* ----------------------------------------------------------------------
   KMC method
   choose and perform an event for site
------------------------------------------------------------------------- */

void AppInbDiff::site_event(int i, class RandomPark *random)
{
  int j,m,n;

  // pick one event from total propensity by accumulating its probability
  // compare prob to threshhold, break when reach it to select event
/*  int rr11 = random->uniform()*5;
  std::cout << "rr11";
  std::cout << rr11 << " ";
  std::cout << "\n";*/
 
  double threshhold = random->uniform() * propensity[i2site[i]];
/*  std::cout << "threshold";
  std::cout << threshhold << " ";
  std::cout << "\n";*/
  double proball = 0.0;

  int ievent = firstevent[i];
  while (1) {
    proball += events[ievent].propensity;
    if (proball >= threshhold) break;
    ievent = events[ievent].next;
  }

  // perform the event
  int rstyle = events[ievent].style;
  int which = events[ievent].which;
  j = events[ievent].jpartner;

/*  std::cout << "which";
  std::cout << which << " ";
  std::cout << "\n";*/
/*  std::cout << "main_total_pvav";
  std::cout << total_pvav << " ";
  std::cout << "\n";*/
    
  //total_time_dt = total_time_dt + 1;
/*  std::cout << "stoptime";
  std::cout << stoptime << " ";
  std::cout << "\n";
  std::cout << "st_time";
  std::cout << time << " ";
  std::cout << "\n";
  std::cout << "dt_step";
  std::cout << dt_step << " ";
  std::cout << "\n";
  
  //if(total_time_dt > 0.5){
  std::cout << "kmc_dt";
  std::cout << dt_kmc << " ";
  std::cout << "\n\n";*/
 //   }
/*  if (which == 36){
  std::cout << "i=";
  std::cout << i << " ";
  std::cout << "\n";
  std::cout << "i_corr_i=";
  std::cout << population[2][i] << " ";
  std::cout << "\n";
  std::cout << "i_corr_j=";
  std::cout << population[3][i] << " ";
  std::cout << "\n";
 
  std::cout << "nbr=";
  std::cout << j << " ";
  std::cout << "\n";
  std::cout << "nbr_corr_i=";
  std::cout << population[2][j] << " ";
  std::cout << "\n";
  std::cout << "nbr_corr_j=";
  std::cout << population[3][j] << " ";
  std::cout << "\n";
 
  std::cout << "pop_i=";
  std::cout << population[6][i] << " ";
  std::cout << "\n";
  
  //if(total_time_dt > 0.5){
  std::cout << "pop_nbr=";
  std::cout << population[6][j]<< " ";
  std::cout << "\n\n";
  }*/
 
  for (int ispecies=0; ispecies < MAX_SPECIES; ispecies++)
    population[ispecies][i] += localDeltaPop[which][ispecies];
  if (rstyle == NBR){
    for (int ispecies=0; ispecies < MAX_SPECIES; ispecies++)
      population[ispecies][j] += nbrDeltaPop[which][ispecies];
  }
  rxn_count[which]++;
 
/* if (which == 36){
  std::cout << "after=";
  std::cout << "i=";
  std::cout << i << " ";
  std::cout << "\n";
  std::cout << "nbr=";
  std::cout << j << " ";
  std::cout << "\n";
  std::cout << "pop_i=";
  std::cout << population[6][i] << " ";
  std::cout << "\n";
  
  //if(total_time_dt > 0.5){
  std::cout << "pop_nbr=";
  std::cout << population[6][j]<< " ";
  std::cout << "\n\n\n";
  }*/
 

/*   if(which == 69) {
      total_pvav = total_pvav -1;
  } 
   if(which == 72) {
      total_pvav = total_pvav -1;
  } 
   if(which == 75) {
      total_pvav = total_pvav -1;
  }
  if(which == 78) {
      total_pvav = total_pvav -1;
  } 
   if(which == 83) {
      total_pvav = total_pvav -1;
  } 
   if(which == 86) {
      total_pvav = total_pvav -1;
  }*/
 
/*  int *clusters_nbr = new int[4]; 
  int *clusters_nbr_nbr = new int[4]; 
  for (int ii=0; ii<4;ii++){
      clusters_nbr[ii] = 10000000;
      clusters_nbr_nbr[ii] = 10000000;
  }*/
 
  // compute propensity changes for participating sites and first neighs
  // ignore update of sites with isite < 0
  // use echeck[] to avoid resetting propensity of same site

  int nsites = 0;

  int isite = i2site[i];
  propensity[isite] = site_propensity(i);
  esites[nsites++] = isite;
  echeck[isite] = 1;

  for (n = 0; n < numneigh[i]; n++) {
    m = neighbor[i][n];
    isite = i2site[m];
    if (isite >= 0 && echeck[isite] == 0) {
      propensity[isite] = site_propensity(m);
      esites[nsites++] = isite;
      echeck[isite] = 1;
    }
  }

  if (rstyle == NBR) {
    for (n = 0; n < numneigh[j]; n++) {
      m = neighbor[j][n];
      isite = i2site[m];
      if (isite >= 0 && echeck[isite] == 0) {
	    propensity[isite] = site_propensity(m);
	    esites[nsites++] = isite;
	    echeck[isite] = 1;
      }
    }
  }
 
 
  solve->update(nsites,esites,propensity);

  // clear echeck array

  for (m = 0; m < nsites; m++) echeck[esites[m]] = 0;
 
  //delete [] clusters_nbr;
  //delete [] clusters_nbr_nbr;

 
}
/* ----------------------------------------------------------------------
  Case:  Centripetal force
   choose and perform an event for site
------------------------------------------------------------------------- */


/*void AppInbDiff::cluster_move(int *clusters_nbr, int *clusters_nbr_nbr)
{
  for (int n=0; n < 4; n++) {
    int m;
    int mm = 0;
    int isite = 0;

    if (clusters_nbr[n] != 10000000) {
      for (int nn = 0; nn < numneigh[clusters_nbr[n]]; nn++) {
        m = neighbor[clusters_nbr[n]][nn];
        isite = i2site[m];
        if (isite >= 0 && echeck[isite] == 0) {
	      propensity[isite] = site_propensity(m);
	      esites[mm++] = isite;
	      echeck[isite] = 1;
        } 
      }
    }
    if (clusters_nbr_nbr[n] != 10000000) {  
      for (int nnn = 0; nnn < numneigh[clusters_nbr_nbr[n]]; nnn++) {
        m = neighbor[clusters_nbr_nbr[n]][nnn];
        isite = i2site[m];
        if (isite >= 0 && echeck[isite] == 0) {
	      propensity[isite] = site_propensity(m);
	      esites[mm++] = isite;
	      echeck[isite] = 1;
        }
      }
    }

    solve->update(mm,esites,propensity);

    for (m = 0; m < mm; m++) echeck[esites[m]] = 0;

  }


  // clear echeck array

}*/

/* ----------------------------------------------------------------------
   clear all events out of list for site I
   add cleared events to free list
------------------------------------------------------------------------- */

void AppInbDiff::clear_events(int i)
{
  int next;
  int index = firstevent[i];
  while (index >= 0) {
    next = events[index].next;
    events[index].next = freeevent;
    freeevent = index;
    nevents--;
    index = next;
  }
  firstevent[i] = -1;
}

/* ----------------------------------------------------------------------
   add an event to list for site I
   event = exchange with site J with probability = propensity
------------------------------------------------------------------------- */

void AppInbDiff::add_event(int i, int rstyle, int which, double propensity,
			  int jpartner)
{
  // grow event list and setup free list

  if (nevents == maxevent) {
    maxevent += DELTAEVENT;
    events = 
      (Event *) memory->srealloc(events,maxevent*sizeof(Event),"app:events");
    for (int m = nevents; m < maxevent; m++) events[m].next = m+1;
    freeevent = nevents;
  }

  int next = events[freeevent].next;

  events[freeevent].style = rstyle;
  events[freeevent].which = which;
  events[freeevent].jpartner = jpartner;
  events[freeevent].propensity = propensity;

  events[freeevent].next = firstevent[i];
  firstevent[i] = freeevent;
  freeevent = next;
  nevents++;
}

/* ---------------------------------------------------------------------- */

void AppInbDiff::add_rxn(int narg, char **arg)
{
  if (narg < 3) error->all(FLERR,"Illegal reaction command");

  // store ID

  if (find_reaction(arg[0]) >= 0) {
    char *str = new char[128];
    sprintf(str,"Reaction ID %s already exists",arg[0]);
    error->all(FLERR,str);
  }

  int n = nreactions + 1;
  rname = (char **) memory->srealloc(rname,n*sizeof(char *),
					  "inb/diff:rname");
  int nlen = strlen(arg[0]) + 1;
  rname[nreactions] = new char[nlen];
  strcpy(rname[nreactions],arg[0]);

  // grow reaction arrays
  memory->grow(localReactants,n,MAX_SPECIES,"inb/diff:localReactants");
  memory->grow(nbrReactants,n,MAX_SPECIES,"inb/diff:nbrReactants");
  memory->grow(localDeltaPop,n,MAX_SPECIES,"inb/diff:localDeltaPop");
  memory->grow(nbrDeltaPop,n,MAX_SPECIES,"inb/diff:nbrDeltaPop");
  memory->grow(rate,n,"inb/diff:rate");
  memory->grow(rxnStyle,n,"inb/diff:rxnStyle");

  for (int ispecies = 0; ispecies < MAX_SPECIES; ispecies++){
    localReactants[nreactions][ispecies] = 0;
    nbrReactants[nreactions][ispecies] = 0;
    localDeltaPop[nreactions][ispecies] = 0;
    nbrDeltaPop[nreactions][ispecies] = 0;
  }

  // find which arg is numeric reaction rate
  char c;
  int rateArg = 1;
  while (rateArg < narg) {
    c = arg[rateArg][0];
    if ((c >= '0' && c <= '9') || c == '+' || c == '-' || c == '.') break;
    rateArg++;
  }

  // find which args are the local and nbr reactants and products
  int localReactArg, nbrReactArg, localProdArg, nbrProdArg;
  localReactArg = nbrReactArg = localProdArg = nbrProdArg = 0;
  for (int iarg = 1; iarg < narg; iarg++) {
    if (strcmp(arg[iarg],"local") == 0){
      if (iarg<rateArg) localReactArg = iarg;
      else localProdArg = iarg;
    } else if (strcmp(arg[iarg],"nbr") == 0){
      if (iarg<rateArg) nbrReactArg = iarg;
      else nbrProdArg = iarg;
    }
  }

  // error checks
  if (rateArg == narg) error->all(FLERR,"Reaction has no numeric rate");
  if (!localReactArg || !nbrReactArg || !localProdArg || !nbrProdArg)
    error->all(FLERR,"One or more 'local' or 'nbr' keywords missing");
  if (localReactArg > nbrReactArg || nbrReactArg > rateArg || rateArg > localProdArg ||
      localProdArg > nbrProdArg)
    error->all(FLERR,"One or more 'local' or 'nbr' keywords out of order");

  // extract reactant and product species names
  int nLocalReactant = 0;
  for (int i = localReactArg+1; i < nbrReactArg; i++) {
    int ispecies = find_species(arg[i]);
    if (ispecies == -1) error->all(FLERR,"Unknown species in reaction command");
    localReactants[nreactions][ispecies]++;
    localDeltaPop[nreactions][ispecies]--;
    nLocalReactant++;
  }

  int nNbrReactant = 0;
  for (int i = nbrReactArg+1; i < rateArg; i++) {
    int ispecies = find_species(arg[i]);
    if (ispecies == -1) error->all(FLERR,"Unknown species in reaction command");
    nbrReactants[nreactions][ispecies]++;
    nbrDeltaPop[nreactions][ispecies]--;
    nNbrReactant++;
  }

  rate[nreactions] = atof(arg[rateArg]);

  int nLocalProduct = 0;
  for (int i = localProdArg+1; i < nbrProdArg; i++) {
    int ispecies = find_species(arg[i]);
    if (ispecies == -1) error->all(FLERR,"Unknown species in reaction command");
    localDeltaPop[nreactions][ispecies]++;
    nLocalProduct++;
  }
  
  int nNbrProduct = 0;
  for (int i = nbrProdArg+1; i < narg; i++) {
    int ispecies = find_species(arg[i]);
    if (ispecies == -1) error->all(FLERR,"Unknown species in reaction command");
    nbrDeltaPop[nreactions][ispecies]++;
    nNbrProduct++;
  }

  // Set the reaction style as local only or not
  rxnStyle[nreactions] = NBR;
  if(!nNbrReactant && !nNbrProduct){
    rxnStyle[nreactions] = LOCAL;
  }

  // additional error checking
  if(!nLocalReactant && !nLocalProduct && !nNbrReactant && !nNbrProduct)
    error->all(FLERR,"Reaction must have at least one reactant or product");
  if(!nLocalReactant && !nLocalProduct)
    error->all(FLERR,"Reaction cannot only act on neighbor");

  nreactions++;
}

/* ---------------------------------------------------------------------- */

void AppInbDiff::add_species(int narg, char **arg)
{
  if (narg == 0) error->all(FLERR,"Illegal species command");

  // grow species arrays

  int n = nspecies + narg;
  sname = (char **) memory->srealloc(sname,n*sizeof(char *),
					  "inb/diff:sname");

  for (int iarg = 0; iarg < narg; iarg++) {
    if (find_species(arg[iarg]) >= 0) {
      char *str = new char[128];
      sprintf(str,"Species ID %s already exists",arg[iarg]);
      error->all(FLERR,str);
    }
    int nlen = strlen(arg[iarg]) + 1;
    sname[nspecies+iarg] = new char[nlen];
    strcpy(sname[nspecies+iarg],arg[iarg]);
  }
  nspecies += narg;
}

/* ----------------------------------------------------------------------
   return reaction index (0 to N-1) for a reaction ID
   return -1 if doesn't exist
------------------------------------------------------------------------- */

int AppInbDiff::find_reaction(char *str)
{
  for (int i = 0; i < nreactions; i++)
    if (strcmp(str,rname[i]) == 0) return i;
  return -1;
}

/* ----------------------------------------------------------------------
   return species index (0 to N-1) for a species ID
   return -1 if doesn't exist
------------------------------------------------------------------------- */

int AppInbDiff::find_species(char *str)
{
  for (int i = 0; i < nspecies; i++)
    if (strcmp(str,sname[i]) == 0) return i;
  return -1;
}

/* ----------------------------------------------------------------------
 * Stub for custom rate multiplier function to be replaced in rxn/diff/custom
------------------------------------------------------------------------- */
double AppInbDiff::custom_multiplier(int i, int rstyle, int which, int jpartner)
{
  return 1.0;
}

