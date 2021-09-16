#ifndef NEUTRINOMFLR_H
#define NEUTRINOMFLR_H

#include "../data/dtypes.h"
#include "../main/main.h"
#include "gadgetconfig.h"

struct nulinear
{
public:
    const double N_nu_massive = 3.0; 
    const double N_nu_eff = 3.0;
    static const int N_tau = 50; // sets the number of neutrino fluids 

    // homogeneous evolution functions
    double tau_t_eV(int t);
    double v2_t_eta_REL(int t, double eta);
    double Ec_t_eta_REL(int t, double eta);
    double eta_convert(double a);
    double Ec_de_eta(double eta);
    double Hc2_Hc02_eta(double eta);
    static double Hubble_Extern(double a);
    double m_nu_eV_parser(void);

    // SuperEasy linear response functions
    double kFs_gen(double k, double a, double m);
    double K_i(double k, double a, double m);
    double poisson_mod_fac(double k, double a);
    double poisson_mod_fac_gen(double k, double a);

//private:
    const double Hc0h = 3.33564095198152e-04;
#define Hc0h2 (Hc0h*Hc0h)
    const double aeta_in = 1e-3;
#define eta_stop (-log(aeta_in))
    
    const double T_CMB_0_K = 2.7255;
   
#define m_nu_eV (93.141613*All.OmegaNu*All.HubbleParam*All.HubbleParam/N_nu_massive)
#define Omega_nu_t_0 (All.OmegaNu/N_tau)
   
#define T_CMB_0_K_4 (T_CMB_0_K*T_CMB_0_K*T_CMB_0_K*T_CMB_0_K)
#define T_NUREL_0_K (0.713765855503608*T_CMB_0_K)
#define m_T_nu (m_nu_eV * 11604.51812 / T_NUREL_0_K)
#define Omega_gam_0 ((4.46911743913795e-07)*T_CMB_0_K_4/(All.HubbleParam*All.HubbleParam))
#define Omega_nurel_0 (0.227107317660239*(N_nu_eff-N_nu_massive)*Omega_gam_0)
#define Omega_rel_0 (Omega_gam_0+Omega_nurel_0)
#define Omega_de_0 (1.0-All.Omega0-All.OmegaNu-Omega_rel_0)
    
    const double w_eos_cdm = 0;
    const double w_eos_gam = 0.33333333333333333;
    const double w0_eos_de = -1.0;
    const double wa_eos_de = 0.0;
    const double wb_eos_de = 0.0;
};

extern nulinear Nulinear;

#endif
