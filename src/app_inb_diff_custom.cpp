/* -------------------------i---------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPPARKS directory.
------------------------------------------------------------------------- */

#include <iostream>
#include <algorithm>

#include "math.h"
#include "string.h"
#include "app_inb_diff_custom.h"
#include "random_park.h"
#include "error.h"

using namespace SPPARKS_NS;

enum{LOCAL,NBR};

/* ---------------------------------------------------------------------- */

AppInbDiffCustom::AppInbDiffCustom(SPPARKS *spk, int narg, char **arg) : 
  AppInbDiff(spk,narg,arg)
{
  // parse arguments for InbDiffCustom class only, not children
  if (strcmp(style,"inb/diff/custom") != 0) return;
  if (narg != 1) error->all(FLERR,"Illegal app_style command");
  /*total_cluster = atof(arg[1]);
  beta = atof(arg[2]);
  radius = atoi(arg[3]);
  ran_com = atof(arg[4]);
  periphery_rate = atof(arg[5]);*/
 
}

/* ----------------------------------------------------------------------
 * Custom rate multiplier function
------------------------------------------------------------------------- */
/*void AppInbDiffCustom::cluster_movement(int i, int rstyle, int which, int jpartner)
{
  if (which==22) {
    int pp;
    for (int ppp=0; ppp<numneigh[i]; ppp++) {
      pp = neighbor[i][ppp];
      if (pp != jpartner){
        int diff = jpartner - i;
        for (int ispecies=0; ispecies < MAX_SPECIES; ispecies++)
          population[ispecies][pp] += localDeltaPop[which][ispecies];
      
        for (int ispecies=0; ispecies < MAX_SPECIES; ispecies++)
          population[ispecies][pp+diff] += nbrDeltaPop[which][ispecies];
      }   
    }
  }
}*/

double AppInbDiffCustom::custom_multiplier(int i, int rstyle, int which, int jpartner)
{
   //if (which==0) { // Diffusion: Activating ligands
 //    if(population[27][i] != 0) {
 //      int bar_fill_ahead = population[18][i]+ population[5][i]+population[6][i]+population[7][i]+population[10][i]+population[11][i]+population[12][i]+population[13][i]+population[19][i] +population[20][i]+population[22][i]+population[23][i]+population[24][i]+population[25][i];
;
  //     if (bar_fill_ahead >=2500) {
   //      return 0.0;
  //     } else return 1.0;
  //   } else return 0.0;
   if (which==2) { // Diffusion: Activating ligands
     if(population[1][jpartner] != 0 && i != jpartner) {
      // int ar_ligand_fill_ahead = population[0][jpartner]+ population[17][jpartner]+population[6][jpartner]+population[7][jpartner]+population[10][jpartner]+population[11][jpartner]+population[12][jpartner]+population[13][jpartner]+population[19][jpartner] +population[20][jpartner]+population[22][jpartner]+population[23][jpartner]+population[24][jpartner]+population[25][jpartner]+population[29][jpartner]+population[30][jpartner]+population[31][jpartner]+population[33][jpartner]+population[34][jpartner];
       //if (ar_ligand_fill_ahead >=2500) {
/*       std::cout << "jpartner=";
       std::cout << jpartner << " ";
       std::cout << "\n";
       std::cout << "i=";
       std::cout << i << " ";
       std::cout << "\n\n";
       std::cout << "pop_j=";
       std::cout << population[1][jpartner] << " ";
       std::cout << "\n";
       std::cout << "pop_i=";
       std::cout << population[1][i]  << " ";
       std::cout << "\n\n";*/
 
         return 1.0;
       } else return 0.0;
    // } else return 0.0;
   } else if (which==3) { // Difusion: Activating receptors
     if(population[1][jpartner] != 0 && i != jpartner ) {
        if(population[15][i] != 0 && population[15][jpartner] != 0 ) {
      // int u_act_fill_ahead = population[5][jpartner]+population[6][jpartner]+population[7][jpartner]+population[10][jpartner]+population[11][jpartner]+population[12][jpartner]+population[13][jpartner]+population[18][jpartner]+population[19][jpartner] +population[20][jpartner]+population[22][jpartner]+population[23][jpartner]+population[24][jpartner]+population[25][jpartner]+population[29][jpartner]+population[30][jpartner]+population[31][jpartner]+population[8][jpartner]+population[33][jpartner]+population[34][jpartner];
      // if (u_act_fill_ahead >=2500) {
            return 1.0;
       } else if (population[15][i] == 0 && population[15][jpartner] == 0 ) {
         return 1.0;
       } else if (population[15][i] != 0 && population[15][jpartner] == 0 ){
         return 0.0;
       } else if (population[15][i] == 0 && population[15][jpartner] != 0 ) {
         return 1.0;
       }
    } else return 0.0;
   } else if (which==4) { //Difussion: Src Kinase
     if(population[1][jpartner] != 0 && i != jpartner) {
      // int ir_inb = population[18][i]+population[19][i] +population[20][i]+population[22][i]+population[23][i]+population[24][i]+population[25][i]+population[33][i]+population[34][i];
      // int u_sfk_fill_ahead = population[8][jpartner]+population[5][jpartner]+population[6][jpartner]+population[7][jpartner]+population[10][jpartner]+population[11][jpartner]+population[12][jpartner]+population[13][jpartner]+population[18][jpartner]+population[19][jpartner] +population[20][jpartner]+population[22][jpartner]+population[23][jpartner]+population[24][jpartner]+population[25][jpartner]+population[29][jpartner]+population[30][jpartner]+population[31][jpartner]+population[33][jpartner]+population[34][jpartner];
      // if (ir_inb <= bir_conc) {
        // if (u_sfk_fill_ahead >=2500) {
           return 1.0;
         } else {
            return 0.0;
        }
   } else if (which==23) { // Diffusion: unbound Inhibitory Receptors
     if(population[16][jpartner] != 0 && i != jpartner) {
         if(population[16][i] != 0 && population[16][jpartner] != 0 ) {
      // int u_act_fill_ahead = population[5][jpartner]+population[6][jpartner]+population[7][jpartner]+population[10][jpartner]+population[11][jpartner]+population[12][jpartner]+population[13][jpartner]+population[18][jpartner]+population[19][jpartner] +population[20][jpartner]+population[22][jpartner]+population[23][jpartner]+population[24][jpartner]+population[25][jpartner]+population[29][jpartner]+population[30][jpartner]+population[31][jpartner]+population[8][jpartner]+population[33][jpartner]+population[34][jpartner];
      // if (u_act_fill_ahead >=2500) {
            return 1.0;
       } else if (population[16][i] == 0 && population[16][jpartner] == 0 ) {
         return 1.0;
       } else if (population[16][i] != 0 && population[16][jpartner] == 0 ){
         return 0.0;
       } else if (population[16][i] == 0 && population[16][jpartner] != 0 ) {
         return 1.0;
       }

      // int micro_inb = population[18][i]+population[19][i] +population[20][i]+population[22][i]+population[23][i]+population[24][i]+population[25][i];
      // int micro_inb_fill_ahead = population[5][jpartner]+population[6][jpartner]+population[7][jpartner]+population[10][jpartner]+population[11][jpartner]+population[12][jpartner]+population[13][jpartner]+population[18][jpartner]+population[19][jpartner] +population[20][jpartner]+population[22][jpartner]+population[23][jpartner]+population[24][jpartner]+population[25][jpartner]+population[29][jpartner]+population[30][jpartner]+population[31][jpartner]+population[8][jpartner]+population[33][jpartner]+population[34][jpartner];
       //if (micro_inb <= bir_conc) {
       //  if( micro_inb_fill_ahead >=2500) {
      } else return 0.0;
      // } else return 0.0;
   //  } else return 0.0;
  } else if (which==44) { // Diffusion: Inhibitory ligands
     if(population[1][jpartner] != 0 && i != jpartner) {
      // int ir_ligand_fill_ahead = population[0][jpartner]+population[17][jpartner] + population[6][jpartner]+population[7][jpartner]+population[10][jpartner]+population[11][jpartner]+population[12][jpartner]+population[13][jpartner]+population[19][jpartner] +population[20][jpartner]+population[22][jpartner]+population[23][jpartner]+population[24][jpartner]+population[25][jpartner]+population[29][jpartner]+population[30][jpartner]+population[31][jpartner]+population[33][jpartner]+population[34][jpartner];

  //     if (ir_ligand_fill_ahead >=2500) {
         return 1.0;
     } else return 0.0;
  } else {

  return 1.0;
// }
}
}

/* ----------------------------------------------------------------------
   print stats
------------------------------------------------------------------------- */

/*void AppInbDiffCustom::stats(int *ir_cluster_center)
{
  delete [] ir_cluster_center;
  char big[8],format[64];
  strcpy(big,BIGINT_FORMAT);

 // FIX: massive hack: food_left updated whenever stats() is called
  int rxn_count_b2B, rxn_count_p2P;
  MPI_Allreduce(&rxn_count[17],&rxn_count_b2B,1,MPI_INT,MPI_SUM,world);
  MPI_Allreduce(&rxn_count[18],&rxn_count_p2P,1,MPI_INT,MPI_SUM,world);
  food_left = num_nutrient - rxn_count_b2B - rxn_count_p2P;

  bigint naccept_all;
  MPI_Allreduce(&naccept,&naccept_all,1,MPI_SPK_BIGINT,MPI_SUM,world);

  sprintf(format,"%%10g %%10%s %%10d %%10d");
  sprintf(strtmp,format,time,naccept_all,nsweeps);
  sprintf(format,"%%10g %%10%s %%10d %%10d",&big[1]);
  sprintf(strtmp,format,time,naccept_all,nsweeps,food_left);

}*/

/* ----------------------------------------------------------------------
   print stats header
------------------------------------------------------------------------- */

/*void AppInbDiffCustom::stats_header(char *strtmp)
{
  //sprintf(strtmp,"%10s %10s %10s %10s","Time","Naccept","Nsweeps","Nfood");
  sprintf(strtmp,"%10s %10s %10s","Time","Naccept","Nsweeps");
}*/
