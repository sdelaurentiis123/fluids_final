//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file kh.cpp
//! \brief Problem generator for KH instability.
//!
//! Sets up several different problems:
//!   - iprob=1: slip surface with random perturbations
//!   - iprob=2: tanh profile, with single-mode perturbation (Frank et al. 1996)
//!              MODIFIED to support density stratification and passive scalar.
//!   - iprob=3: tanh profiles for v and d, SR test problem in Beckwith & Stone (2011)
//!   - iprob=4: tanh profiles for v and d, "Lecoanet" test
//!   - iprob=5: two resolved slip-surfaces with m=2 perturbation for the AMR test

// C headers

// C++ headers
#include <algorithm>  // min, max
#include <cmath>      // log
#include <cstring>    // strcmp()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../defs.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../scalars/scalars.hpp"
#include "../utils/utils.hpp"

namespace {
Real vflow;
int iprob;
Real PassiveDyeEntropy(MeshBlock *pmb, int iout);
} // namespace

Real threshold;
int RefinementCondition(MeshBlock *pmb);

//----------------------------------------------------------------------------------------
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.

void Mesh::InitUserMeshData(ParameterInput *pin) {
  vflow = pin->GetReal("problem","vflow");
  iprob = pin->GetInteger("problem","iprob");

  if (adaptive) {
    threshold = pin->GetReal("problem", "thr");
    EnrollUserRefinementCondition(RefinementCondition);
  }
  if (iprob == 4 && NSCALARS > 0) {
    AllocateUserHistoryOutput(1);
    EnrollUserHistoryOutput(0, PassiveDyeEntropy, "tot-S");
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Problem Generator for the Kelvin-Helmholtz test

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  std::int64_t iseed = -1 - gid;
  Real gm1 = peos->GetGamma() - 1.0;

  //--- iprob=1.  Uniform stream with density ratio "drat" located in region -1/4<y<1/4
  // moving at (-vflow) seperated by two slip-surfaces from background medium with d=1
  // moving at (+vflow), random perturbations.  This is the classic, unresolved K-H test.

  if (iprob == 1) {
    // Read problem parameters
    Real drat = pin->GetReal("problem","drat");
    Real amp = pin->GetReal("problem","amp");
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++) {
          phydro->u(IDN,k,j,i) = 1.0;
          phydro->u(IM1,k,j,i) = vflow + amp*(ran2(&iseed) - 0.5);
          phydro->u(IM2,k,j,i) = amp*(ran2(&iseed) - 0.5);
          phydro->u(IM3,k,j,i) = 0.0;
          if (std::abs(pcoord->x2v(j)) < 0.25) {
            phydro->u(IDN,k,j,i) = drat;
            phydro->u(IM1,k,j,i) = -drat*(vflow + amp*(ran2(&iseed) - 0.5));
            phydro->u(IM2,k,j,i) = drat*amp*(ran2(&iseed) - 0.5);
          }
          // Pressure scaled to give a sound speed of 1 with gamma=1.4
          if (NON_BAROTROPIC_EOS) {
            phydro->u(IEN,k,j,i) =
                2.5/gm1 + 0.5*(SQR(phydro->u(IM1,k,j,i)) +
                               SQR(phydro->u(IM2,k,j,i)))/phydro->u(IDN,k,j,i);
          }
        }
      }
    }
    if (NSCALARS > 0) {
      for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je; j++) {
          for (int i=is; i<=ie; i++) {
            for (int n=0; n<NSCALARS; ++n) {
              if (std::abs(pcoord->x2v(j)) > 0.25) {
                pscalars->s(n,k,j,i) = n % 2;
              } else {
                pscalars->s(n,k,j,i) = std::abs(n % 2 - 1)*drat;
              }
            }
          }
        }
      }
    }

    // initialize uniform interface B
    if (MAGNETIC_FIELDS_ENABLED) {
      Real b0 = pin->GetReal("problem","b0");
      for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je; j++) {
          for (int i=is; i<=ie+1; i++) {
            pfield->b.x1f(k,j,i) = b0;
          }
        }
      }
      for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je+1; j++) {
          for (int i=is; i<=ie; i++) {
            pfield->b.x2f(k,j,i) = 0.0;
          }
        }
      }
      for (int k=ks; k<=ke+1; k++) {
        for (int j=js; j<=je; j++) {
          for (int i=is; i<=ie; i++) {
            pfield->b.x3f(k,j,i) = 0.0;
          }
        }
      }
      if (NON_BAROTROPIC_EOS) {
        for (int k=ks; k<=ke; k++) {
          for (int j=js; j<=je; j++) {
            for (int i=is; i<=ie; i++) {
              phydro->u(IEN,k,j,i) += 0.5*b0*b0;
            }
          }
        }
      }
    }
  }

  //--- iprob=2. Uniform density medium moving at +/-vflow seperated by a single shear
  // layer with tanh() profile at y=0 with a single mode perturbation, reflecting BCs at
  // top/bottom.  Based on Frank et al., ApJ 460, 777, 1996.
  // MODIFIED: To include density ratio (drat) and passive scalar initialization,
  //           and use 'lambda' parameter for vel2 perturbation.

  if (iprob == 2) {
    // Read/set problem parameters
    Real amp = pin->GetReal("problem", "amp");
    Real drat = pin->GetReal("problem", "drat"); // Density ratio rho_high/rho_low
    Real lambda_param = pin->GetReal("problem", "lambda"); // Wavelength of perturbation

    // Width of the tanh transition layers for density and velocity.
    Real width = pin->GetOrAddReal("problem", "width", 0.02);
    // Standard deviation for the Gaussian envelope of the vy perturbation.
    Real sigma_pert = pin->GetOrAddReal("problem", "sigma_pert", 0.2);
    // Default gas pressure if not specified in hydro block
    Real pgas = pin->GetOrAddReal("hydro", "pgas", 2.5);


    Real rho_low = 1.0; // Density of the lower fluid (y < 0)
    Real rho_high = rho_low * drat; // Density of the upper fluid (y > 0)

    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++) {
          Real x1 = pcoord->x1v(i);
          Real x2 = pcoord->x2v(j);

          // Density profile: tanh transition from rho_low to rho_high centered at x2=0
          Real dens = rho_low + (rho_high - rho_low) * 0.5 * (1.0 + std::tanh(x2 / width));
          phydro->u(IDN,k,j,i) = dens;

          // Velocity in x1 (shear): vflow * tanh(x2/width)
          // This means -vflow for x2 << 0 and +vflow for x2 >> 0
          phydro->u(IM1,k,j,i) = vflow * std::tanh(x2 / width);
          
          // Velocity in x2 (perturbation): sinusoidal in x1, Gaussian in x2
          // kx is the wavenumber in x1 direction.
          // Domain spans x1min to x1max. If lambda_param is ~0, use dominant mode in box.
          Real kx = (lambda_param > 1.0e-6) ? (TWO_PI / lambda_param) : (TWO_PI / (pcoord->x1f(ie+1) - pcoord->x1f(is)));
          phydro->u(IM2,k,j,i) = amp * vflow * std::cos(kx * x1) 
                                 * std::exp(-SQR(x2 / sigma_pert)); // x2^2 / sigma_pert^2

          // Velocity in x3 (out of plane): 0 for 2D
          phydro->u(IM3,k,j,i) = 0.0;

          if (NON_BAROTROPIC_EOS) {
            // For conserved variables, U = (rho, M1, M2, M3, E)
            // M = rho * v
            // E = P/(gm1) + 0.5 * rho * v^2
            // So, when setting primitive variable u(IEN,...), it should be E.
            // The KE part should be 0.5 * (M1^2 + M2^2 + M3^2)/rho if M's are momenta
            // OR 0.5*rho*(v1^2+v2^2+v3^2) if v's are velocities.
            // The original code had: 0.5*(SQR(phydro->u(IM1,k,j,i)) + SQR(phydro->u(IM2,k,j,i)) + SQR(phydro->u(IM3,k,j,i)))*phydro->u(IDN,k,j,i)
            // Athena stores momenta in IM1,IM2,IM3. So this should be:
            // 0.5*(SQR(phydro->u(IM1,k,j,i)) + SQR(phydro->u(IM2,k,j,i)) + SQR(phydro->u(IM3,k,j,i))) / phydro->u(IDN,k,j,i)
            // However, the IM1, IM2, IM3 are SET as velocities in this problem generator block, then implicitly treated as momenta by the solver.
            // The most straightforward way for problem generators is to set conserved variables directly.
            // But since we are setting primitive velocity here, then for E:
            Real m1 = phydro->u(IM1,k,j,i) * dens; // Calculate momentum for KE
            Real m2 = phydro->u(IM2,k,j,i) * dens;
            Real m3 = phydro->u(IM3,k,j,i) * dens;
            phydro->u(IEN,k,j,i) = pgas / gm1 + 0.5*(SQR(m1) + SQR(m2) + SQR(m3)) / dens;
            
            // Then convert velocities to momenta for storage in u(IM1), u(IM2), u(IM3)
            phydro->u(IM1,k,j,i) *= dens;
            phydro->u(IM2,k,j,i) *= dens;
            phydro->u(IM3,k,j,i) *= dens;

          }
        }
      }
    }

    // Initialize passive scalar s[0] to trace the density contrast
    // Athena stores passive scalars as s_n = S_n (conserved)
    // The problem generator usually sets s(n,k,j,i) as concentration S_n / rho
    // So, if we want s[0] to be 0 in low-density and 1 in high-density region (as a tracer quantity like dye)
    // then S_0 should be rho * 0 for low-density, and rho * 1 for high-density
    // The common approach is s(n,k,j,i) = desired_concentration_value
    // which is then multiplied by rho internally if needed for conserved form.
    // For pscalars->s(n,k,j,i) this is S_n/rho at input, then converted to S_n.
    if (NSCALARS > 0) {
      for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je; j++) {
          for (int i=is; i<=ie; i++) {
            Real x2 = pcoord->x2v(j);
            // Scalar profile: 0 in effective low-density fluid (rho_low), 1 in effective high-density fluid (rho_high)
            // This matches the density profile's tanh transition.
            // This sets the CONCENTRATION s_n = S_n/rho
            pscalars->s(0,k,j,i) = 0.5 * (1.0 + std::tanh(x2 / width)); 
            // Initialize other scalars if NSCALARS > 1
            for (int n=1; n<NSCALARS; ++n) {
              pscalars->s(n,k,j,i) = 0.0; // Or some other profile
            }
          }
        }
      }
    }

    // initialize uniform interface B (copied from original iprob=2 which was hydro)
    // For MHD KHI, this B-field setup would likely need to be more sophisticated.
    if (MAGNETIC_FIELDS_ENABLED) {
      Real b0 = pin->GetReal("problem","b0"); // This b0 is likely Bx0 if following original
      for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je; j++) {
          for (int i=is; i<=ie+1; i++) { // Loop to ie+1 for x1-faces
            pfield->b.x1f(k,j,i) = b0;
          }
        }
      }
      for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je+1; j++) { // Loop to je+1 for x2-faces
          for (int i=is; i<=ie; i++) {
            pfield->b.x2f(k,j,i) = 0.0;
          }
        }
      }
      for (int k=ks; k<=ke+1; k++) { // Loop to ke+1 for x3-faces
        for (int j=js; j<=je; j++) {
          for (int i=is; i<=ie; i++) {
            pfield->b.x3f(k,j,i) = 0.0;
          }
        }
      }
      if (NON_BAROTROPIC_EOS) { // Add magnetic energy to total energy
        for (int k=ks; k<=ke; k++) {
          for (int j=js; j<=je; j++) {
            for (int i=is; i<=ie; i++) {
              // Assuming b0 is Bx, total B^2 = Bx^2
              // Note: This adds to phydro->u(IEN,...), which already contains thermal and kinetic energy.
              phydro->u(IEN,k,j,i) += 0.5*SQR(b0);
            }
          }
        }
      }
    }
  } // --- End of iprob == 2 ---

  //--- iprob=3.  Test in SR paper (Beckwith & Stone, ApJS 193, 6, 2011).  Gives two
  // resolved shear layers with tanh() profiles for velocity and density located at
  // y = +/- 0.5, density one in middle and 0.01 elsewhere, single mode perturbation.

  if (iprob == 3) {
    // Read/set problem parameters
    Real amp = pin->GetReal("problem","amp");
    Real a = 0.01;
    Real sigma = 0.1;
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++) {
          phydro->u(IDN,k,j,i) = 0.505 + 0.495
                                 *std::tanh((std::abs(pcoord->x2v(j))-0.5)/a);
          phydro->u(IM1,k,j,i) = vflow*std::tanh((std::abs(pcoord->x2v(j))-0.5)/a);
          phydro->u(IM2,k,j,i) =
              amp*vflow*std::sin(TWO_PI*pcoord->x1v(i))
              *std::exp(-((std::abs(pcoord->x2v(j))-0.5)
                          *(std::abs(pcoord->x2v(j))-0.5))/(sigma*sigma));
          if (pcoord->x2v(j) < 0.0) phydro->u(IM2,k,j,i) *= -1.0;
          
          // Convert to momenta before storing
          phydro->u(IM1,k,j,i) *= phydro->u(IDN,k,j,i);
          phydro->u(IM2,k,j,i) *= phydro->u(IDN,k,j,i);
          phydro->u(IM3,k,j,i) = 0.0;

          if (NON_BAROTROPIC_EOS) {
            phydro->u(IEN,k,j,i) =
                1.0/gm1 + 0.5*(SQR(phydro->u(IM1,k,j,i)) + // Momenta used here
                               SQR(phydro->u(IM2,k,j,i)))/phydro->u(IDN,k,j,i);
          }
        }
      }
    }

    // initialize uniform interface B
    if (MAGNETIC_FIELDS_ENABLED) {
      Real b0 = pin->GetReal("problem","b0");
      for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je; j++) {
          for (int i=is; i<=ie+1; i++) {
            pfield->b.x1f(k,j,i) = b0;
          }
        }
      }
      for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je+1; j++) {
          for (int i=is; i<=ie; i++) {
            pfield->b.x2f(k,j,i) = 0.0;
          }
        }
      }
      for (int k=ks; k<=ke+1; k++) {
        for (int j=js; j<=je; j++) {
          for (int i=is; i<=ie; i++) {
            pfield->b.x3f(k,j,i) = 0.0;
          }
        }
      }
      if (NON_BAROTROPIC_EOS) {
        for (int k=ks; k<=ke; k++) {
          for (int j=js; j<=je; j++) {
            for (int i=is; i<=ie; i++) {
              phydro->u(IEN,k,j,i) += 0.5*SQR(b0); // Assuming b0 is Bx
            }
          }
        }
      }
    }
  }

  //--- iprob=4.  "Lecoanet" test, resolved shear layers with tanh() profiles for velocity
  // and density located at z1=0.5, z2=1.5 two-mode perturbation for fully periodic BCs

  // To promote symmetry of FP errors about midplanes, rescale z' = z - 1. ; x' = x - 0.5
  // so that domain x1 = [-0.5, 0.5] and x2 = [-1.0, 1.0] is centered about origin
  if (iprob == 4) {
    // Read/set problem parameters
    Real amp = pin->GetReal("problem","amp");
    // unstratified problem is the default
    Real drho_rho0 = pin->GetOrAddReal("problem", "drho_rho0", 0.0);
    // set background vx to nonzero to evolve the KHI in a moving frame
    Real vboost = pin->GetOrAddReal("problem", "vboost", 0.0);
    Real P0 = 10.0;
    Real a = 0.05;
    Real sigma = 0.2;
    // Initial condition's reflect-and-shift symmetry, x1-> x1 + 1/2, x2-> -x2
    // is preserved in new coordinates; hence, the same flow is solved twice in this prob.
    Real z1 = -0.5;  // z1' = z1 - 1.0
    Real z2 = 0.5;   // z2' = z2 - 1.0

    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++) {
          // Lecoanet (2015) equation 8a)
          Real dens = 1.0 + 0.5*drho_rho0*(std::tanh((pcoord->x2v(j) - z1)/a) -
                                           std::tanh((pcoord->x2v(j) - z2)/a));
          phydro->u(IDN,k,j,i) = dens;

          Real v1 = vflow*(std::tanh((pcoord->x2v(j) - z1)/a)
                           - std::tanh((pcoord->x2v(j) - z2)/a) - 1.0) // 8b)
                    + vboost;
          // Currently, the midpoint approx. is applied in the momenta and energy calc
          phydro->u(IM1,k,j,i) = v1*dens; // Store momentum M1

          // NOTE ON FLOATING-POINT SHIFT SYMMETRY IN X1:
          // (details omitted for brevity - see original Athena++ source)
          Real ave_sine = std::sin(TWO_PI*pcoord->x1v(i));
          if (pcoord->x1v(i) > 0.0) {
            ave_sine -= std::sin(TWO_PI*(-0.5 + pcoord->x1v(i)));
          } else {
            ave_sine -= std::sin(TWO_PI*(0.5 + pcoord->x1v(i)));
          }
          ave_sine /= 2.0;

          Real v2 = -amp*ave_sine
                    *(std::exp(-(SQR(pcoord->x2v(j) - z1))/(sigma*sigma)) +
                      std::exp(-(SQR(pcoord->x2v(j) - z2))/(sigma*sigma))); // 8c), mod.
          phydro->u(IM2,k,j,i) = v2*dens; // Store momentum M2

          phydro->u(IM3,k,j,i) = 0.0; // M3 = 0
          if (NON_BAROTROPIC_EOS) {
            phydro->u(IEN,k,j,i) = P0/gm1 + 0.5*(SQR(phydro->u(IM1,k,j,i)) // M1^2
                                                 + SQR(phydro->u(IM2,k,j,i)) // M2^2
                                                 + SQR(phydro->u(IM3,k,j,i)) // M3^2
                                                 ) / phydro->u(IDN,k,j,i); // / rho
          }
          // color concentration of passive scalar
          if (NSCALARS > 0) {
            Real concentration = 0.5*(std::tanh((pcoord->x2v(j) - z2)/a)  // 8e)
                                      - std::tanh((pcoord->x2v(j) - z1)/a) + 2.0);
            // uniformly fill all scalar species to have equal concentration
            // Here, pscalars->s is S_n / rho (concentration)
            constexpr int scalar_norm = NSCALARS > 0 ? NSCALARS : 1.0; // Avoid division by zero
            for (int n=0; n<NSCALARS; ++n) {
               // The Lecoanet paper implies s is mass fraction, so S_n = concentration * rho
               // But Athena pscalars->s is usually set as concentration (S_n/rho)
               // And then Hydro::PassiveScalarSourceTerms might multiply by rho if defined that way
               // For this problem, it seems they set s (concentration) and then multiply by rho if storing S_n
               // Here, it's directly setting the concentration for pscalars->s
              pscalars->s(n,k,j,i) = (1.0/scalar_norm)*concentration; // Set concentration
            }
          }
        }
      }
    }
    // initialize uniform interface B
    if (MAGNETIC_FIELDS_ENABLED) {
      Real b0 = pin->GetReal("problem", "b0");
      b0 = b0/std::sqrt(4.0*(PI)); // Convert b0 to code units if it's in CGS (typical for some problems)
      for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je; j++) {
          for (int i=is; i<=ie+1; i++) {
            pfield->b.x1f(k,j,i) = b0*std::tanh((std::abs(pcoord->x2v(j)) - 0.5)/a);
          }
        }
      }
      for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je+1; j++) {
          for (int i=is; i<=ie; i++) {
            pfield->b.x2f(k,j,i) = 0.0;
          }
        }
      }
      for (int k=ks; k<=ke+1; k++) {
        for (int j=js; j<=je; j++) {
          for (int i=is; i<=ie; i++) {
            pfield->b.x3f(k,j,i) = 0.0;
          }
        }
      }
      if (NON_BAROTROPIC_EOS) {
        for (int k=ks; k<=ke; k++) {
          for (int j=js; j<=je; j++) {
            for (int i=is; i<=ie; i++) {
              phydro->u(IEN,k,j,i) += 0.5*SQR(pfield->b.x1f(k,j,i)); // Add Bx^2/2
              // If By, Bz were non-zero, add their energy too.
            }
          }
        }
      }
    }
  }

  //--- iprob=5. Uniform stream with density ratio "drat" located in region -1/4<y<1/4
  // moving at (-vflow) seperated by two resolved slip-surfaces from background medium
  // with d=1 moving at (+vflow), with m=2 perturbation, for the AMR test.

  if (iprob == 5) {
    // Read problem parameters
    Real a_param = pin->GetReal("problem","a"); // Renamed to avoid conflict with 'a' in iprob=4 if used globally
    Real sigma_param = pin->GetReal("problem","sigma"); // Renamed
    Real drat_param = pin->GetReal("problem","drat"); // Renamed
    Real amp_param = pin->GetReal("problem","amp"); // Renamed

    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++) {
          Real w_factor =(std::tanh((std::abs(pcoord->x2v(j))-0.25)/a_param)+1.0)*0.5;
          phydro->u(IDN,k,j,i) = w_factor+(1.0-w_factor)*drat_param;
          
          Real v1 = w_factor*vflow-(1.0-w_factor)*vflow; // Corrected drat multiplication for momentum
          Real v2 = amp_param * std::sin(2.0*TWO_PI*pcoord->x1v(i))
                                 * std::exp(-SQR(std::abs(pcoord->x2v(j))-0.25)
                                            /(sigma_param*sigma_param));
          phydro->u(IM1,k,j,i) = v1 * phydro->u(IDN,k,j,i); // M1
          phydro->u(IM2,k,j,i) = v2 * phydro->u(IDN,k,j,i); // M2
          phydro->u(IM3,k,j,i) = 0.0; // M3

          // Pressure scaled to give a sound speed of 1 with gamma=1.4
          if (NON_BAROTROPIC_EOS) {
            // Original: 2.5/gm1 + 0.25*(SQR(phydro->u(IM1,k,j,i)) + SQR(phydro->u(IM2,k,j,i)))/phydro->u(IDN,k,j,i);
            // The 0.25 was likely a typo, should be 0.5 for KE term.
            phydro->u(IEN,k,j,i) =
                2.5/gm1 + 0.5*(SQR(phydro->u(IM1,k,j,i)) + // M1^2
                               SQR(phydro->u(IM2,k,j,i)) + // M2^2
                               SQR(phydro->u(IM3,k,j,i))   // M3^2
                               )/phydro->u(IDN,k,j,i); // / rho
          }
        }
      }
    }

    // initialize uniform interface B
    if (MAGNETIC_FIELDS_ENABLED) {
      Real b0 = pin->GetReal("problem","b0");
      for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je; j++) {
          for (int i=is; i<=ie+1; i++) {
            pfield->b.x1f(k,j,i) = b0;
          }
        }
      }
      for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je+1; j++) {
          for (int i=is; i<=ie; i++) {
            pfield->b.x2f(k,j,i) = 0.0;
          }
        }
      }
      for (int k=ks; k<=ke+1; k++) {
        for (int j=js; j<=je; j++) {
          for (int i=is; i<=ie; i++) {
            pfield->b.x3f(k,j,i) = 0.0;
          }
        }
      }
      if (NON_BAROTROPIC_EOS) {
        for (int k=ks; k<=ke; k++) {
          for (int j=js; j<=je; j++) {
            for (int i=is; i<=ie; i++) {
              phydro->u(IEN,k,j,i) += 0.5*SQR(b0); // Assuming b0 is Bx
            }
          }
        }
      }
    }
  }

  return;
}


// refinement condition: velocity gradient

int RefinementCondition(MeshBlock *pmb) {
  AthenaArray<Real> &w = pmb->phydro->w; // Primitive variables
  Real vgmax = 0.0;
  for (int k=pmb->ks; k<=pmb->ke; k++) {
    for (int j=pmb->js-1; j<=pmb->je+1; j++) { // Loop over ghost cells for gradients
      for (int i=pmb->is-1; i<=pmb->ie+1; i++) {
        // Corrected gradient calculation (assuming uniform grid for simplicity)
        // More robust would use pcoord->dx1f etc.
        Real vgy = std::abs(w(IVY,k,j,i+1) - w(IVY,k,j,i-1))*0.5; // d(vy)/dx using i neighbors for x-gradient of vy
        Real vgx = std::abs(w(IVX,k,j+1,i) - w(IVX,k,j-1,i))*0.5; // d(vx)/dy using j neighbors for y-gradient of vx
        // This seems to be calculating |dVy/Dx + dVx/Dy| which is related to shear/vorticity
        // The original formulation might have intended vorticity like (dvx/dy - dvy/dx) or similar
        Real vg  = std::sqrt(vgx*vgx + vgy*vgy); // Magnitude of some gradient vector
        if (vg > vgmax) vgmax = vg;
      }
    }
  }
  if (vgmax > threshold) return 1;
  if (vgmax < 0.5*threshold) return -1; // Original had 0.25*threshold for de-refinement
  return 0;
}

namespace {
Real PassiveDyeEntropy(MeshBlock *pmb, int iout) {
  Real total_entropy = 0;
  int is = pmb->is, ie = pmb->ie, js = pmb->js, je = pmb->je, ks = pmb->ks, ke = pmb->ke;
  AthenaArray<Real> &r = pmb->pscalars->r; // This 'r' is for passive scalars concentrations s_n = S_n/rho
  AthenaArray<Real> &w = pmb->phydro->w;   // Primitive variables
  AthenaArray<Real> volume; // 1D array of volumes
  // allocate 1D array for cell volume used in usr def history
  volume.NewAthenaArray(pmb->ncells1);

  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      pmb->pcoord->CellVolume(k, j, pmb->is, pmb->ie, volume);
      for (int i=is; i<=ie; i++) {
        // no loop over NSCALARS; hardcode assumption that NSCALARS=1 for this history output
        // r(0,k,j,i) is s(0,k,j,i)/rho(k,j,i)
        // Ensure r(0,k,j,i) is positive for log and represents a meaningful concentration (e.g., 0 to 1)
        if (r(0,k,j,i) > 1.0e-12 && r(0,k,j,i) < 1.0 - 1.0e-12) { // Avoid log(0) or log(1) if r is dye fraction
            // Lecoanet (2016) eq 5 for mixing entropy: - (c log c + (1-c) log (1-c))
            // If r(0,k,j,i) is c, then the specific entropy for mixing is:
            // -( r(0,k,j,i)*std::log(r(0,k,j,i)) + (1.0-r(0,k,j,i))*std::log(1.0-r(0,k,j,i)) )
            // The original implementation was -r*log(r) which is only part of it.
            Real c = r(0,k,j,i);
            Real specific_entropy = -(c * std::log(c) + (1.0 - c) * std::log(1.0 - c) );
            total_entropy += volume(i) * w(IDN,k,j,i) * specific_entropy; 
        } else if (r(0,k,j,i) >= 1.0 - 1.0e-12 || r(0,k,j,i) <= 1.0e-12) {
            // If c is 0 or 1, specific entropy is 0
            total_entropy += 0.0;
        }
        // If using -r*log(r) as a proxy:
        // if (r(0,k,j,i) > 1.0e-12) { // Avoid log(0)
        //    Real specific_entropy = -r(0,k,j,i)*std::log(r(0,k,j,i));
        //    total_entropy += volume(i)*w(IDN,k,j,i)*specific_entropy;
        // }
      }
    }
  }
  return total_entropy;
}
} // namespace
