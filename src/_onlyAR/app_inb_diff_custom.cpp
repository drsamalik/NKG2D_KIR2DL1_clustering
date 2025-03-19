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
  if (narg != 7) error->all(FLERR,"Illegal app_style command");
  total_cluster = atof(arg[1]);
  beta = atof(arg[2]);
  radius = atoi(arg[3]);
  center_vicinity = atoi(arg[4]);
 // frac_vel = atof(arg[5]);
  ran_com = atof(arg[5]);
  ran_com1 = atof(arg[6]);
//  bir_conc = atoi(arg[7]);
 
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
     if(population[26][jpartner] != 0) {
       int ar_ligand_fill_ahead = population[0][jpartner]+ population[17][jpartner]+population[6][jpartner]+population[7][jpartner]+population[10][jpartner]+population[11][jpartner]+population[12][jpartner]+population[13][jpartner]+population[19][jpartner] +population[20][jpartner]+population[22][jpartner]+population[23][jpartner]+population[24][jpartner]+population[25][jpartner]+population[29][jpartner]+population[30][jpartner]+population[31][jpartner]+population[33][jpartner]+population[34][jpartner];
       if (ar_ligand_fill_ahead >=2500) {
         return 0.0;
       } else return 1.0;
     } else return 0.0;
   } else if (which==3) { // Difusion: Activating receptors
     if(population[26][jpartner] != 0) {
       int u_act_fill_ahead = population[5][jpartner]+population[6][jpartner]+population[7][jpartner]+population[10][jpartner]+population[11][jpartner]+population[12][jpartner]+population[13][jpartner]+population[18][jpartner]+population[19][jpartner] +population[20][jpartner]+population[22][jpartner]+population[23][jpartner]+population[24][jpartner]+population[25][jpartner]+population[29][jpartner]+population[30][jpartner]+population[31][jpartner]+population[8][jpartner]+population[33][jpartner]+population[34][jpartner];
       if (u_act_fill_ahead >=2500) {
         return 0.0;
       } else return 1.0;
     } else return 0.0;
   } else if (which==4) { //Difussion: Src Kinase
     if(population[26][jpartner] != 0) {
      // int ir_inb = population[18][i]+population[19][i] +population[20][i]+population[22][i]+population[23][i]+population[24][i]+population[25][i]+population[33][i]+population[34][i];
       int u_sfk_fill_ahead = population[8][jpartner]+population[5][jpartner]+population[6][jpartner]+population[7][jpartner]+population[10][jpartner]+population[11][jpartner]+population[12][jpartner]+population[13][jpartner]+population[18][jpartner]+population[19][jpartner] +population[20][jpartner]+population[22][jpartner]+population[23][jpartner]+population[24][jpartner]+population[25][jpartner]+population[29][jpartner]+population[30][jpartner]+population[31][jpartner]+population[33][jpartner]+population[34][jpartner];
      // if (ir_inb <= bir_conc) {
         if (u_sfk_fill_ahead >=2500) {
           return 0.0;
         } else {
            return 1.0;
        }
      // } else return 0.0; 
     } else return 0.0;
/*   } else if (which==555) { // Diffusion: Vav (cytosol)
     if(population[26][jpartner] != 0) {
       int u_vav_fill_ahead = population[9][jpartner]+population[10][jpartner]+population[11][jpartner]+population[13][jpartner]+population[14][jpartner]+population[24][jpartner]+population[25][jpartner];
       if (u_vav_fill_ahead >=300000) {
         return 0.0;
       } else return 1.0;
     } else return 0.0;*/
  } else if (which==5) { // Activating receptors + SFK
     if(population[1][i] != 0) {
         return 1.0;
     } else return 0.0;
   } else if (which==17) { // Diffusion:free  pVav (cytosol)
     if(population[26][jpartner] != 0) {
       int p_vav_fill_ahead = population[9][jpartner]+population[10][jpartner]+population[11][jpartner]+population[13][jpartner]+population[14][jpartner]+population[24][jpartner]+population[25][jpartner]+population[28][jpartner]+population[29][jpartner]+population[30][jpartner]+population[31][jpartner]+population[32][jpartner]+population[33][jpartner]+population[34][jpartner];
       if (p_vav_fill_ahead >=300000) {
         return 0.0;
       } else return 1.0;
     } else return 0.0;
   } else if (which==56) { // Influx: Activating receptors
     if(population[16][i] != 0 && population[27][jpartner] != 0) {
       int u_act_fill_ahead = population[5][jpartner]+population[6][jpartner]+population[7][jpartner]+population[10][jpartner]+population[11][jpartner]+population[12][jpartner]+population[13][jpartner]+population[18][jpartner]+population[19][jpartner] +population[20][jpartner]+population[22][jpartner]+population[23][jpartner]+population[24][jpartner]+population[25][jpartner]+population[29][jpartner]+population[30][jpartner]+population[31][jpartner]+population[8][jpartner]+population[33][jpartner]+population[34][jpartner];
       if (u_act_fill_ahead >=2500) {
         return 0.0;
       } else return 1.0;
     } else return 0.0;
   } else if (which==57) { // Diffusion: unbound Inhibitory Receptors
     if(population[26][jpartner] != 0) {
       int micro_inb = population[18][i]+population[19][i] +population[20][i]+population[22][i]+population[23][i]+population[24][i]+population[25][i];
       int micro_inb_fill_ahead = population[5][jpartner]+population[6][jpartner]+population[7][jpartner]+population[10][jpartner]+population[11][jpartner]+population[12][jpartner]+population[13][jpartner]+population[18][jpartner]+population[19][jpartner] +population[20][jpartner]+population[22][jpartner]+population[23][jpartner]+population[24][jpartner]+population[25][jpartner]+population[29][jpartner]+population[30][jpartner]+population[31][jpartner]+population[8][jpartner]+population[33][jpartner]+population[34][jpartner];
       //if (micro_inb <= bir_conc) {
         if( micro_inb_fill_ahead >=2500) {
           return 0.0;
         } else return 1.0;
      // } else return 0.0;
     } else return 0.0;
   } else if (which==58) { // Diffusion: Inhibitory ligands
     if(population[26][jpartner] != 0) {
       int ir_ligand_fill_ahead = population[0][jpartner]+population[17][jpartner] + population[6][jpartner]+population[7][jpartner]+population[10][jpartner]+population[11][jpartner]+population[12][jpartner]+population[13][jpartner]+population[19][jpartner] +population[20][jpartner]+population[22][jpartner]+population[23][jpartner]+population[24][jpartner]+population[25][jpartner]+population[29][jpartner]+population[30][jpartner]+population[31][jpartner]+population[33][jpartner]+population[34][jpartner];
;
       if (ir_ligand_fill_ahead >=2500) {
         return 0.0;
       } else return 1.0;
     } else return 0.0;
  } else if (which==61) { //  Inhibitory receptors + SFK
     if(population[1][i] != 0) {
         return 1.0;
     } else return 0.0;
/*   } else if (which==399) { // Diffusion: SHP (cytosol)
     if(population[26][jpartner] != 0) {
       int u_shp_fill_ahead =  population[21][jpartner]+population[22][jpartner]+population[24][jpartner]+population[25][jpartner]+population[30][jpartner];//+population[14][jpartner]+population[23][jpartner]+population[24][jpartner]+population[20][jpartner]+population[21][jpartner];
       if (u_shp_fill_ahead >=700000) {
         return 0.0;
       } else return 1.0;
     } else return 0.0;*/
   } else if (which>=28 && which<=33) { // AR microcluster formation - pVAV dependent
    if ( (population[1][jpartner] != 0) && (population[1][i] != 0 )) {
      double energy_voxel = -((population[11][i]+population[14][i]+population[24][i]+population[25][i]+population[31][i]+population[32][i]+population[33][i]+population[34][i])*(population[11][i]+population[14][i]+population[24][i]+population[25][i]+population[31][i]+population[32][i]+population[33][i]+population[34][i]));
      double energy_voxel_min = 0.0;

      double energy_nbr_voxel = -((population[11][jpartner]+population[14][jpartner]+population[24][jpartner]+population[25][jpartner]+population[31][jpartner]+population[32][jpartner]+population[33][jpartner]+population[34][jpartner])*(population[11][jpartner]+population[14][jpartner]+population[24][jpartner]+population[25][jpartner]+population[31][jpartner]+population[32][jpartner]+population[33][jpartner]+population[34][jpartner]));
/*        std::cout << "\n";
        std::cout << "energy_voxel=";
        std::cout << energy_voxel << " ";
        std::cout << "\n";
        std::cout << "energy_nbr_voxel=";
        std::cout << energy_nbr_voxel << " ";
        std::cout << "\n";
        std::cout << "i_x";
        std::cout << population[2][i] << " ";
        std::cout << "\n";
        std::cout << "i_y";
        std::cout << population[3][i] << " ";
        std::cout << "\n";
        std::cout << "j_x";
        std::cout << population[2][jpartner] << " ";
        std::cout << "\n";
        std::cout << "j_y";
        std::cout << population[3][jpartner] << " ";
        std::cout << "\n";
        std::cout << "i_value";
        std::cout << i << " ";
        std::cout << "\n";
        std::cout << "j_value";
        std::cout << jpartner << " ";
        std::cout << "\n";*/

      double energy_nbr_voxel_min = 0.0;
      int k;
      int jj;
      for (int kk=0; kk < numneigh[i]; kk++) {
        k = neighbor[i][kk];
        double energy_voxel_temp = -((population[11][k]+population[14][k]+population[24][k]+population[25][k]+population[31][k]+population[32][k]+population[33][k]+population[34][k])*(population[11][k]+population[14][k]+population[24][k]+population[25][k]+population[31][k]+population[32][k]+population[33][k]+population[34][k]));
        energy_voxel_min = energy_voxel_temp;
/*        std::cout << "\n";
        std::cout << "k_i_nbr_value";
        std::cout << k << " ";
        std::cout << "\n";
        std::cout << "energy_voxel_temp";
        std::cout << energy_voxel_temp << " ";
        std::cout << "\n";*/
        jj = neighbor[jpartner][kk];
        double energy_nbr_voxel_temp = -((population[11][jj]+population[14][jj]+population[24][jj]+population[25][jj]+population[31][jj]+population[32][jj]+population[33][jj]+population[34][jj])*(population[11][jj]+population[14][jj]+population[24][jj]+population[25][jj]+population[31][jj]+population[32][jj]+population[33][jj]+population[34][jj]));
        energy_nbr_voxel_min = energy_nbr_voxel_temp;
/*        std::cout << "\n";
        std::cout << "jj_j_nbr_value";
        std::cout << jj << " ";
        std::cout << "\n";
        std::cout << "energy_nbr_voxel_temp";
        std::cout << energy_nbr_voxel_temp << " ";*/
 
        if (energy_voxel > energy_voxel_min) {
          energy_voxel = energy_voxel_min;
        } 
        if (energy_nbr_voxel > energy_nbr_voxel_min) {
          energy_nbr_voxel = energy_nbr_voxel_min;
        }  
      }
/*        std::cout << "\n";
        std::cout << "energy_voxel_min";
        std::cout << energy_voxel << " ";
        std::cout << "\n";
        std::cout << "energy_nbr_voxel_min";
        std::cout << energy_nbr_voxel << " ";*/
 
      double vav_quorum;
      double vav_dep_pot;
      double total_vav_dep_pot;
      if (energy_voxel == 0 && energy_nbr_voxel == 0) {
        vav_quorum = 0.0;
      } else if (energy_voxel !=0 || energy_nbr_voxel !=0 ) {
        vav_dep_pot = (beta*(energy_voxel-energy_nbr_voxel));
        total_vav_dep_pot = exp(vav_dep_pot)/(1.0+exp(vav_dep_pot));
        vav_quorum = std::min(1.0,total_vav_dep_pot);
/*        std::cout << "\n";
        std::cout << "vav_dep_pot";
        std::cout << vav_dep_pot << " ";
        std::cout << "\n";
        std::cout << "total_vav_dep_pot";
        std::cout << total_vav_dep_pot << " ";
         std::cout << "\n";
        std::cout << "vav_quorum";
        std::cout << vav_quorum << " ";*/
 
      }
       int pot_u_act_fill_ahead = population[5][jpartner]+population[6][jpartner]+population[7][jpartner]+population[10][jpartner]+population[11][jpartner]+population[12][jpartner]+population[13][jpartner]+population[18][jpartner]+population[19][jpartner] +population[20][jpartner]+population[22][jpartner]+population[23][jpartner]+population[24][jpartner]+population[25][jpartner]+population[29][jpartner]+population[30][jpartner]+population[31][jpartner]+population[8][jpartner]+population[33][jpartner]+population[34][jpartner];
      if (pot_u_act_fill_ahead >=2500) {
        return 0.0;
      } else return vav_quorum;
    } else return 0.0;
   } else if (which>=34 && which<=36) { // AR microcluster formation - pVAV dependent
    if ( (population[1][jpartner] != 0) && (population[1][i] != 0 )) {
      double energy_voxel_vav_prime = -((population[11][i]+population[14][i]+population[24][i]+population[25][i]+population[31][i]+population[32][i]+population[33][i]+population[34][i])*(population[11][i]+population[14][i]+population[24][i]+population[25][i]+population[31][i]+population[32][i]+population[33][i]+population[34][i]));
      double energy_voxel_min_vav_prime = 0.0;

      double energy_nbr_voxel_vav_prime = -((population[11][jpartner]+population[14][jpartner]+population[24][jpartner]+population[25][jpartner]+population[31][jpartner]+population[32][jpartner]+population[33][jpartner]+population[34][jpartner])*(population[11][jpartner]+population[14][jpartner]+population[24][jpartner]+population[25][jpartner]+population[31][jpartner]+population[32][jpartner]+population[33][jpartner]+population[34][jpartner]));

      double energy_nbr_voxel_min_vav_prime = 0.0;
      int k;
      int jj;
      for (int kk=0; kk < numneigh[i]; kk++) {
        k = neighbor[i][kk];
        double energy_voxel_temp_vav_prime = -((population[11][k]+population[14][k]+population[24][k]+population[25][k]+population[31][k]+population[32][k]+population[33][k]+population[34][k])*(population[11][k]+population[14][k]+population[24][k]+population[25][k]+population[31][k]+population[32][k]+population[33][k]+population[34][k]));
        energy_voxel_min_vav_prime = energy_voxel_temp_vav_prime;
        jj = neighbor[jpartner][kk];
        double energy_nbr_voxel_temp_vav_prime = -((population[11][jj]+population[14][jj]+population[24][jj]+population[25][jj]+population[31][jj]+population[32][jj]+population[33][jj]+population[34][jj])*(population[11][jj]+population[14][jj]+population[24][jj]+population[25][jj]+population[31][jj]+population[32][jj]+population[33][jj]+population[34][jj]));
        energy_nbr_voxel_min_vav_prime = energy_nbr_voxel_temp_vav_prime;
        if (energy_voxel_vav_prime > energy_voxel_min_vav_prime) {
          energy_voxel_vav_prime = energy_voxel_min_vav_prime;
        } 
        if (energy_nbr_voxel_vav_prime > energy_nbr_voxel_min_vav_prime) {
          energy_nbr_voxel_vav_prime = energy_nbr_voxel_min_vav_prime;
        }  
      }
      double vav_quorum_vav_prime;
      double vav_dep_pot_vav_prime;
      double total_vav_dep_pot_vav_prime;
      if (energy_voxel_vav_prime == 0 && energy_nbr_voxel_vav_prime == 0) {
        vav_quorum_vav_prime = 0.0;
      } else if (energy_voxel_vav_prime !=0 || energy_nbr_voxel_vav_prime !=0 ) {
        vav_dep_pot_vav_prime = (beta*(energy_voxel_vav_prime - energy_nbr_voxel_vav_prime));
        total_vav_dep_pot_vav_prime = exp(vav_dep_pot_vav_prime)/(1.0+exp(vav_dep_pot_vav_prime));
        vav_quorum_vav_prime = std::min(1.0,total_vav_dep_pot_vav_prime);
      }
       int pot_u_act_fill_ahead_vav_prime = population[28][jpartner]+population[29][jpartner]+population[30][jpartner]+population[31][jpartner]+population[32][jpartner]+population[33][jpartner]+population[34][jpartner];
       int pot_u_act_fill_ahead_vav_prime_ar = population[5][jpartner]+population[6][jpartner]+population[7][jpartner]+population[10][jpartner]+population[11][jpartner]+population[12][jpartner]+population[13][jpartner]+population[18][jpartner]+population[19][jpartner] +population[20][jpartner]+population[22][jpartner]+population[23][jpartner]+population[24][jpartner]+population[25][jpartner]+population[29][jpartner]+population[30][jpartner]+population[31][jpartner]+population[8][jpartner]+population[33][jpartner]+population[34][jpartner];
      if ((pot_u_act_fill_ahead_vav_prime_ar < 2500) && (pot_u_act_fill_ahead_vav_prime < 250)) {
        return vav_quorum_vav_prime;
      } else return 0.0;
    } else return 0.0;
  } else if (which>=37 && which<=42) { //Centripetal movement of AR microclusters
      if ((population[1][jpartner] != 0) && (population[1][i] !=0 )){
        double quorum_ar = 0.0;
        int coor_i;
        int coor_j;
        coor_i = population[2][i];
        coor_j = population[3][i];
        int dia = radius*2;
        int width = (dia/2.0)-center_vicinity;
       
// centrial movement
        if( (coor_i < dia && coor_j < width) || (coor_i < width && coor_j < dia) || (coor_i > (dia-width) && coor_j < dia ) || (coor_i < dia && coor_j > dia-width) ) {
          if ( jpartner == i-1) {
           quorum_ar = ran_com* ((1+(((total_pvav)/(total_cluster+total_pvav))*((coor_i/((radius-2.0)*1.0))-1)))/4.0) +(1-ran_com)*0.25;
          } else if ( jpartner == i+1) {
           quorum_ar = ran_com*( (1-(((total_pvav)/(total_cluster+total_pvav))*((coor_i/((radius-2.0)*1.0))-1)))/4.0) + (1-ran_com)*0.25;
          } else if ((jpartner - i) < -1) {
           quorum_ar = ran_com*( (1+(((total_pvav)/(total_cluster+total_pvav))*((coor_j/((radius-2.0)*1.0))-1)))/4.0) +(1-ran_com)*0.25;
          } else if ((jpartner - i) > 1) {
           quorum_ar = ran_com*( (1-(((total_pvav)/(total_cluster+total_pvav))*((coor_j/((radius-2.0)*1.0))-1)))/4.0) +(1-ran_com)*0.25;
          }
        } else {
            if ( jpartner == i-1) {
            quorum_ar = ran_com1*( (1+((((total_pvav)/((total_cluster)+total_pvav)))*((coor_i/((radius-2.0)*1.0))-1)))/4.0) + (1-ran_com1)*0.25;
          } else if ( jpartner == i+1) {
            quorum_ar = ran_com1*( (1-((((total_pvav)/((total_cluster)+total_pvav)))*((coor_i/((radius-2.0)*1.0))-1)))/4.0) + (1-ran_com1)*0.25;
          } else if ((jpartner - i) < -1) {
            quorum_ar = ran_com1*( (1+((((total_pvav)/((total_cluster)+total_pvav)))*((coor_j/((radius-2.0)*1.0))-1)))/4.0) + (1-ran_com1)*0.25;
          } else if ((jpartner - i) > 1) {
            quorum_ar = ran_com1*( (1-((((total_pvav)/((total_cluster)+total_pvav)))*((coor_j/((radius-2.0)*1.0))-1)))/4.0) + (1-ran_com1)*0.25;
          }
        }
        int nbr;
        int diff = jpartner-i;
        
        int micro_act_fill_ahead = population[5][jpartner]+population[6][jpartner]+population[7][jpartner]+population[10][jpartner]+population[11][jpartner]+population[12][jpartner]+population[13][jpartner]+population[18][jpartner]+population[19][jpartner] +population[20][jpartner]+population[22][jpartner]+population[23][jpartner]+population[24][jpartner]+population[25][jpartner]+population[29][jpartner]+population[30][jpartner]+population[31][jpartner]+population[8][jpartner]+population[33][jpartner]+population[34][jpartner];
        
        int *nbr_i_fill = new int[4];
        for (int kk =0; kk<4;kk++) {
            nbr_i_fill[kk] = 1000000;
        }
        int lnbr_i_fill = 0;
        int rnbr_i_fill = 0;
        int dnbr_i_fill = 0;
        int unbr_i_fill = 0;
        if (which == 37) { 
            if (population[1][i] != 0) {
                lnbr_i_fill =population[6][i];
            }
        } else if (which == 38) {  
            if (population[1][i] != 0) {
                lnbr_i_fill =population[12][i];
            }
        } else if (which == 39) {  
            if (population[1][i] != 0) {
                lnbr_i_fill =population[7][i];
            }
        } else if (which == 40) {  
            if (population[1][i] != 0) {
                lnbr_i_fill =population[10][i];
            }
        } else if (which == 41) {  
            if (population[1][i] != 0) {
                lnbr_i_fill =population[13][i];
            }
        } else if (which == 42) {  
            if (population[1][i] != 0) {
                lnbr_i_fill =population[11][i];
            }
        }
 
// check the AR vacancy for first 4 nearest nbr voxels 
        if(lnbr_i_fill >= 1) {
            nbr_i_fill[0] = i;
        }
 
       //if ((micro_act_fill_ahead < 2500) && (nbr_i_fill[0] == i)) {
       if ((micro_act_fill_ahead >= 2500)) {
                delete [] nbr_i_fill;
                  return 0.0;
        } else {
                  delete [] nbr_i_fill;
                return quorum_ar;
        }
    } else return 0.0;
  } else if (which >= 43 && which <=45) { //Centripetal movement of AR attached to vav_prime microclusters
      if ((population[1][jpartner] != 0) && (population[1][i] !=0 )){
        double quorum_vav_prime_ar = 0.0;
        int coor_ii;
        int coor_jj;
        coor_ii = population[2][i];
        coor_jj = population[3][i];
        int dia_vav_prime_ar = radius*2;
        int width_vav_prime_ar = (dia_vav_prime_ar/2.0)-center_vicinity;
       
// centrial movement
        if( (coor_ii < dia_vav_prime_ar && coor_jj < width_vav_prime_ar) || (coor_ii < width_vav_prime_ar && coor_jj < dia_vav_prime_ar) || (coor_ii > (dia_vav_prime_ar - width_vav_prime_ar) && coor_jj < dia_vav_prime_ar ) || (coor_ii < dia_vav_prime_ar && coor_jj > dia_vav_prime_ar - width_vav_prime_ar) ) {
          if ( jpartner == i-1) {
           quorum_vav_prime_ar = ran_com* ((1+(((total_pvav)/(total_cluster+total_pvav))*((coor_ii/((radius-2.0)*1.0))-1)))/4.0) +(1-ran_com)*0.25;
          } else if ( jpartner == i+1) {
           quorum_vav_prime_ar = ran_com*( (1-(((total_pvav)/(total_cluster+total_pvav))*((coor_ii/((radius-2.0)*1.0))-1)))/4.0) + (1-ran_com)*0.25;
          } else if ((jpartner - i) < -1) {
           quorum_vav_prime_ar = ran_com*( (1+(((total_pvav)/(total_cluster+total_pvav))*((coor_jj/((radius-2.0)*1.0))-1)))/4.0) +(1-ran_com)*0.25;
          } else if ((jpartner - i) > 1) {
           quorum_vav_prime_ar = ran_com*( (1-(((total_pvav)/(total_cluster+total_pvav))*((coor_jj/((radius-2.0)*1.0))-1)))/4.0) +(1-ran_com)*0.25;
          }
        } else {
            if ( jpartner == i-1) {
            quorum_vav_prime_ar = ran_com1*( (1+((((total_pvav)/((total_cluster)+total_pvav)))*((coor_ii/((radius-2.0)*1.0))-1)))/4.0) + (1-ran_com1)*0.25;
          } else if ( jpartner == i+1) {
            quorum_vav_prime_ar = ran_com1*( (1-((((total_pvav)/((total_cluster)+total_pvav)))*((coor_ii/((radius-2.0)*1.0))-1)))/4.0) + (1-ran_com1)*0.25;
          } else if ((jpartner - i) < -1) {
            quorum_vav_prime_ar = ran_com1*( (1+((((total_pvav)/((total_cluster)+total_pvav)))*((coor_jj/((radius-2.0)*1.0))-1)))/4.0) + (1-ran_com1)*0.25;
          } else if ((jpartner - i) > 1) {
            quorum_vav_prime_ar = ran_com1*( (1-((((total_pvav)/((total_cluster)+total_pvav)))*((coor_jj/((radius-2.0)*1.0))-1)))/4.0) + (1-ran_com1)*0.25;
          }
        }
        int nbr;
        int diff_vav_prime_ar = jpartner-i;
        
        int micro_vav_prime_ar_fill_ahead = population[28][jpartner]+population[29][jpartner]+population[30][jpartner]+population[31][jpartner]+population[32][jpartner]+population[33][jpartner]+population[34][jpartner];
        
        int micro_vav_prime_ar_ar_fill_ahead = population[5][jpartner]+population[6][jpartner]+population[7][jpartner]+population[10][jpartner]+population[11][jpartner]+population[12][jpartner]+population[13][jpartner]+population[18][jpartner]+population[19][jpartner] +population[20][jpartner]+population[22][jpartner]+population[23][jpartner]+population[24][jpartner]+population[25][jpartner]+population[29][jpartner]+population[30][jpartner]+population[31][jpartner]+population[8][jpartner]+population[33][jpartner]+population[34][jpartner];
        int *nbr_i_fill_vav_prime_ar = new int[4];
        for (int kk =0; kk<4;kk++) {
            nbr_i_fill_vav_prime_ar[kk] = 1000000;
        }
        int lnbr_i_fill_vav_prime_ar = 0;
        int rnbr_i_fill_vav_prime_ar = 0;
        int dnbr_i_fill_vav_prime_ar = 0;
        int unbr_i_fill_vav_prime_ar = 0;
        if (which == 43) {
            if (population[1][i] != 0) {
                lnbr_i_fill_vav_prime_ar =population[29][i];
            }
        } else if (which == 44) {
            if (population[1][i] != 0) {
                lnbr_i_fill_vav_prime_ar =population[30][i];
            }
        } else if (which == 45) {
            if (population[1][i] != 0) {
                lnbr_i_fill_vav_prime_ar =population[31][i];
            }
        }

 
// check the AR vacancy for first 4 nearest nbr voxels 
        if(lnbr_i_fill_vav_prime_ar >= 1) {
            nbr_i_fill_vav_prime_ar[0] = i;
        }

        //if ((micro_vav_prime_ar_fill_ahead < 250) && (micro_vav_prime_ar_ar_fill_ahead < 2500) && (nbr_i_fill_vav_prime_ar[0] == i)) {
        if ((micro_vav_prime_ar_fill_ahead < 250) && (micro_vav_prime_ar_ar_fill_ahead < 2500)) {
                delete [] nbr_i_fill_vav_prime_ar;
                return quorum_vav_prime_ar;
         } else {
                 delete [] nbr_i_fill_vav_prime_ar;
                 return 0.0;
         }
    } else return 0.0;
  } else if (which >= 46 && which <=47) { //Centripetal movement of vav_prime microclusters
      if ((population[1][jpartner] != 0) && (population[1][i] !=0 )){
        double quorum_vav_prime = 0.0;
        int coor_ii;
        int coor_jj;
        coor_ii = population[2][i];
        coor_jj = population[3][i];
        int dia_vav_prime = radius*2;
        int width_vav_prime = (dia_vav_prime/2.0)-center_vicinity;
       
// centrial movement
        if( (coor_ii < dia_vav_prime && coor_jj < width_vav_prime) || (coor_ii < width_vav_prime && coor_jj < dia_vav_prime) || (coor_ii > (dia_vav_prime - width_vav_prime) && coor_jj < dia_vav_prime ) || (coor_ii < dia_vav_prime && coor_jj > dia_vav_prime - width_vav_prime) ) {
          if ( jpartner == i-1) {
           quorum_vav_prime = ran_com* ((1+(((total_pvav)/(total_cluster+total_pvav))*((coor_ii/((radius-2.0)*1.0))-1)))/4.0) +(1-ran_com)*0.25;
          } else if ( jpartner == i+1) {
           quorum_vav_prime = ran_com*( (1-(((total_pvav)/(total_cluster+total_pvav))*((coor_ii/((radius-2.0)*1.0))-1)))/4.0) + (1-ran_com)*0.25;
          } else if ((jpartner - i) < -1) {
           quorum_vav_prime = ran_com*( (1+(((total_pvav)/(total_cluster+total_pvav))*((coor_jj/((radius-2.0)*1.0))-1)))/4.0) +(1-ran_com)*0.25;
          } else if ((jpartner - i) > 1) {
           quorum_vav_prime = ran_com*( (1-(((total_pvav)/(total_cluster+total_pvav))*((coor_jj/((radius-2.0)*1.0))-1)))/4.0) +(1-ran_com)*0.25;
          }
        } else {
            if ( jpartner == i-1) {
            quorum_vav_prime = ran_com1*( (1+((((total_pvav)/((total_cluster)+total_pvav)))*((coor_ii/((radius-2.0)*1.0))-1)))/4.0) + (1-ran_com1)*0.25;
          } else if ( jpartner == i+1) {
            quorum_vav_prime = ran_com1*( (1-((((total_pvav)/((total_cluster)+total_pvav)))*((coor_ii/((radius-2.0)*1.0))-1)))/4.0) + (1-ran_com1)*0.25;
          } else if ((jpartner - i) < -1) {
            quorum_vav_prime = ran_com1*( (1+((((total_pvav)/((total_cluster)+total_pvav)))*((coor_jj/((radius-2.0)*1.0))-1)))/4.0) + (1-ran_com1)*0.25;
          } else if ((jpartner - i) > 1) {
            quorum_vav_prime = ran_com1*( (1-((((total_pvav)/((total_cluster)+total_pvav)))*((coor_jj/((radius-2.0)*1.0))-1)))/4.0) + (1-ran_com1)*0.25;
          }
        }
        int nbr;
        int diff_vav_prime = jpartner-i;
        
        int micro_vav_prime_fill_ahead = population[28][jpartner]+population[29][jpartner]+population[30][jpartner]+population[31][jpartner]+population[32][jpartner]+population[33][jpartner]+population[34][jpartner];
        
        int *nbr_i_fill_vav_prime = new int[4];
        for (int kk =0; kk<4;kk++) {
            nbr_i_fill_vav_prime[kk] = 1000000;
        }
        int lnbr_i_fill_vav_prime = 0;
        int rnbr_i_fill_vav_prime = 0;
        int dnbr_i_fill_vav_prime = 0;
        int unbr_i_fill_vav_prime = 0;
        if (which == 46) {
            if (population[1][i] != 0) {
                lnbr_i_fill_vav_prime =population[28][i];
            }
        } else if (which == 47) {
            if (population[1][i] != 0) {
                lnbr_i_fill_vav_prime =population[32][i];
            }
        }

// check the AR vacancy for first 4 nearest nbr voxels 
        if(lnbr_i_fill_vav_prime >= 1) {
            nbr_i_fill_vav_prime[0] = i;
        }

        //if ((micro_vav_prime_fill_ahead < 250) && (nbr_i_fill_vav_prime[0] == i)) {
        if ((micro_vav_prime_fill_ahead >= 250)) {
                delete [] nbr_i_fill_vav_prime;
                return 0.0;
         } else {
                 delete [] nbr_i_fill_vav_prime;
                return quorum_vav_prime;
         }
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
