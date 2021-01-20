/////////////////////////////////////////////////////////////////////////////
//
// FILE NAME
//   teapotbase.hh
//
// AUTHOR
//   Jeff Holmes, ORNL, jzh@ornl.gov
//   Joshua Abrams, Knox College, jabrams@knox.edu
//   Steven Bunch, University of Tennessee, sbunch2@utk.edu
//
// Modified by Andrei Shishlo
//   12/30/05
//   03/27/15 useCharge added. By default the charge from Bunch will be used.
//            if useCharge != 1 the charge will be assumed +1
//
// Checked by Jeff Holmes
//   02/2012
//
// DESCRIPTION
//   Define elementary functions for different elements
//
/////////////////////////////////////////////////////////////////////////////
#ifndef TEAPOT_BASE_H
#define TEAPOT_BASE_H

#include "Bunch.hh"
#include "OU_Function.hh"

namespace teapot_base
{
    void init_factorial();
    void delete_factorial();

    void rotatexy(Bunch* bunch, double anglexy);

    void drifti(Bunch* bunch, int i, double length);
    void drift(Bunch* bunch, double length);

    void multpi(Bunch* bunch, int i, int pole, double kl, int skew, int useCharge);
    void multp(Bunch* bunch, int pole, double kl, int skew, int useCharge);

    void multpfringeIN(Bunch* bunch, int pole, double kl, int skew, int useCharge);
    void multpfringeOUT(Bunch* bunch, int pole, double kl, const int skew, int useCharge);

	void wrapbunch(Bunch* bunch, double length);

    void kick(Bunch* bunch, double kx, double ky, double kE, int useCharge);

    void quad1(Bunch* bunch, double length, double kq, int useCharge);
    void quad2(Bunch* bunch, double length);
    void quad3(Bunch* bunch, double length, double kq, int useCharge);

    void quadfringeIN(Bunch* bunch, double kq, int useCharge);
    void quadfringeOUT(Bunch* bunch, double kq, int useCharge);

    void wedgerotate(Bunch* bunch, double e, int frinout);
    void wedgedrift(Bunch* bunch, double e, int inout);
    void wedgebend(Bunch* bunch, double e, int inout, double rho, int nsteps);

    void bend1(Bunch* bunch, double length, double th);
    void bend2(Bunch* bunch, double length);
    void bend3(Bunch* bunch, double th);
    void bend4(Bunch* bunch, double th);

    void bendfringeIN(Bunch* bunch, double rho);
    void bendfringeOUT(Bunch* bunch, double rho);

    void soln(Bunch* bunch, double length, double B, int useCharge);

    void wedgebendCF(Bunch* bunch, double e, int inout,
                     double rho,
                     int vecnum,
                     std::vector<int>& pole,
                     std::vector<double>& kl,
                     std::vector<int>& skew,
                     int nsteps, int useCharge);

    void RingRF(Bunch* bunch, double ring_length, int harmonic_numb, double voltage, double phase_s, int useCharge);
    
    void dipoleGeneralKick(Bunch* bunch, double effLength, double strength,double fieldDirection);
    void dipoleXKick(Bunch* bunch, double effLength, double strength);
    void dipoleYKick(Bunch* bunch, double effLength, double strength);
    void dipoleGeneralKickStrip(Bunch* bunch, Bunch* failedToStripBunch, OrbitUtils::Function* survivalProbFunction, OrbitUtils::Function* inverseFunction, double effLength, double strength,double fieldDirection);
    
}

#endif  //TEAPOT_BASE_H

