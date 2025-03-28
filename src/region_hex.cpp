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

#include "stdlib.h"
#include "string.h"
#include "region_hex.h"
#include "domain.h"
#include "error.h"

using namespace SPPARKS_NS;

#define BIG 1.0e20

/* ---------------------------------------------------------------------- */

RegHex::RegHex(SPPARKS *spk, int narg, char **arg) : Region(spk, narg, arg)
{
  options(narg-8,&arg[8]);
  xy = 0;
  if (strcmp(arg[2],"INF") == 0 || strcmp(arg[2],"EDGE") == 0) {
    if (domain->box_exist == 0) 
      error->all(FLERR,"Cannot use region INF or EDGE when box does not exist");
    if (strcmp(arg[2],"INF") == 0) xlo = -BIG;
    else xlo = domain->boxxlo;
  } else xlo = xscale*atof(arg[2]);

  if (strcmp(arg[3],"INF") == 0 || strcmp(arg[3],"EDGE") == 0) {
    if (domain->box_exist == 0) 
      error->all(FLERR,"Cannot use region INF or EDGE when box does not exist");
    if (strcmp(arg[3],"INF") == 0) xhi = BIG;
    else xhi = domain->boxxhi;
  } else xhi = xscale*atof(arg[3]);

  if (strcmp(arg[4],"INF") == 0 || strcmp(arg[4],"EDGE") == 0) {
      error->all(FLERR,"Cannot use region INF or EDGE in y direction with style HEX");
  } else {
    ylo = yscale*atof(arg[4]);
    xy -= atof(arg[4]);
  }

  if (strcmp(arg[5],"INF") == 0 || strcmp(arg[5],"EDGE") == 0) {
      error->all(FLERR,"Cannot use region INF or EDGE in y direction with style HEX");
  } else {
    yhi = yscale*atof(arg[5]);
    xy += atof(arg[5]);
  }

  if (strcmp(arg[6],"INF") == 0 || strcmp(arg[6],"EDGE") == 0) {
    if (domain->box_exist == 0) 
      error->all(FLERR,"Cannot use region INF or EDGE when box does not exist");
    if (strcmp(arg[6],"INF") == 0) zlo = -BIG;
    else zlo = domain->boxzlo;
  } else zlo = zscale*atof(arg[6]);

  if (strcmp(arg[7],"INF") == 0 || strcmp(arg[7],"EDGE") == 0) {
    if (domain->box_exist == 0) 
      error->all(FLERR,"Cannot use region INF or EDGE when box does not exist");
    if (strcmp(arg[7],"INF") == 0) zhi = BIG;
    else zhi = domain->boxzhi;
  } else zhi = zscale*atof(arg[7]);

  xy = xyscale*xy;

  // error check

  if (xlo > xhi || ylo > yhi || zlo > zhi || xy < 0)
    error->all(FLERR,"Illegal region block command. All hi must be greater than paired lo; xy >=0.");

  // extent of block
  
  extent_xlo = xlo;
  extent_xhi = xhi;
  extent_ylo = ylo;
  extent_yhi = yhi;
  extent_zlo = zlo;
  extent_zhi = zhi;
  extent_xy = xy; //by definition this is a distance relative to xlo
}

/* ---------------------------------------------------------------------- */

int RegHex::match(double x, double y, double z)
{
  int inside;
  double prdx = xhi-xlo; 
  if (z >= zlo && z <= zhi && x >= xlo && x <= (xhi+xy) && y >= ylo && y <= yhi) {
    double delx = x - xlo;
    double dely = y - ylo;
    double reffrac = (yhi-ylo)/xy;
    double Rfrac = dely/(delx-prdx);
    double Lfrac = dely/(delx);
    if ((delx >= xy && (delx <= (prdx))) || ((delx <= xy) && (Lfrac <= reffrac)) || ((delx > (prdx)) && (Rfrac >= reffrac))) inside = 1;
    else inside = 0;

  } else inside = 0;
  return !(inside ^ interior);         // 1 if same, 0 if different
}
