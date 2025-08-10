///////////////////////////////////////////////////////////////////////////////                        
/// \class OscProb::PMNS_OQS                                                                                 
///                                                                                                        
/// \brief Implements neutrino oscillations using an open quantum system approach.                         
///                                                                                                        
/// Implementation of three-flavor neutrino oscillations in vacuum and matter,                        
/// incorporating quantum dissipation effects.                                                          
/// Based on the formalism described in: https://doi.org/10.48550/arXiv.hep-ph/0208166.               
///                                                                                              
/// This developement is part of the QGRANT project (ID: 101068013),                                    
/// founded by the HORIZON-MSCA-2021-PF-01-01 programme.                                              
///                                                                                               
/// \author Alba Domi - alba.domi\@fau.de                                                                
/// \author Joao Coelho - jcoelho\@apc.in2p3.fr                                                     
///////////////////////////////////////////////////////////////////////////////


#ifndef PMNS_OQS_H
#define PMNS_OQS_H

#include <Eigen/Core>

#include "PMNS_DensityMatrix.h"

namespace OscProb {

  class PMNS_OQS : public PMNS_DensityMatrix {
    public:
      PMNS_OQS();          ///< Constructor
      virtual ~PMNS_OQS(); ///< Destructor

    static constexpr int SU3_DIM = 9;  ///< Dimension of the SU(3) Gell-Mann basis.

    virtual void SetPower(int n);  ///< Set power-law index for the energy dependence of decoherence parameters.
    virtual void SetDecoElement(int i, double val);  ///< Set value of the a_i decoherence element in Gell-Mann basis.
    virtual void SetDecoAngle(int i, int j, double th);  ///< Set mixing angle between two decoherence parameters a_i, a_j.

    virtual int    GetPower();  ///< Get the currently set power-law index for decoherence parameters.
    virtual double GetDecoElement(int i);  ///< Get the value of the a_i decoherence parameter in Gell-Mann basis.
    virtual double GetDecoAngle(int i, int j);  ///< Get mixing angle between two decoherence parameters \
a_i, a_j.                                               

    virtual double GetHGM(int i, int j);  ///< Get element i,j of the dissipator matrix in Gell-Mann basis.
    virtual double GetDissipatorElement(int i, int j);  ///< Get dissipator element in Gell-Mann basis.

      virtual void    SetIsNuBar(bool isNuBar);
      virtual matrixD ProbMatrix(int nflvi, int nflvf);

    protected:
    virtual void BuildHms();  ///< Build standard Hamiltonian in  mass basis (Hms).
    virtual void SetHeff(NuPath p);  ///< Set Effective Hamiltonian in Vacuum Mass Basis (VMB).
    virtual void SetHGM();  ///< Set Effective Hamiltonian in Gell-Mann Basis.
    virtual void SetDissipator();  ///< Set the dissipator in Gell-Mann Basis.
    virtual void SetM();  ///< Build the matrix M used in evolution equations.
    virtual void RotateState(bool to_mass); ///< Rotate rho to/from mass basis.
    virtual void ChangeBaseToGM();  ///< Change base to Gell-Mann.
    virtual void ChangeBaseToSU3();  ///< Change base to SU(3).

    virtual void Propagate();
    /// Propagate neutrino through a single path
    virtual void PropagatePath(NuPath p);

    int     fPower;  ///< Power-law index (n).
      vectorD fa;   ///< a vector
      matrixD fcos; ///< cosines ai . aj

    vectorD fR;  ///< Initial state of density matrix in Gell-Mann basis.
    vectorD fRt;  ///< Time evolution of density matrix in Gell-Mann basis. 
    matrixC fUM; ///< PMNS Matrix
    matrixC fHeff;  ///< Effective Hamiltonian in VMB (3x3).
    matrixD fHGM;  ///< Effective Hamiltonian in GMB (9x9).
    matrixD fD; ///< Off-diagonal, 9x9 dissipator

      Eigen::MatrixXd fM; ///< Buffer matrix for exponential

      bool fBuiltDissipator; ///< Flag to rebuilt D

  };

} // namespace OscProb

#endif

///////////////////////////////////////////////////////////////////////////////
