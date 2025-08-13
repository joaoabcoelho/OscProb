///////////////////////////////////////////////////////////////////////////////
/// \class OscProb::PMNS_OQS
///
/// \brief Implements neutrino oscillations using an open quantum system
/// approach.
///
/// Implementation of three-flavor neutrino oscillations in vacuum and matter,
/// incorporating quantum dissipation effects.
/// Based on the formalism described in:
/// https://doi.org/10.48550/arXiv.hep-ph/0208166.
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

      /// Dimension of the SU(3) Gell-Mann basis.
      static constexpr int SU3_DIM = 9;

      /// Set power-law index for the energy dependence of decoherence
      /// parameters.
      virtual void SetPower(int n);
      /// Set value of the a_i decoherence element in Gell-Mann basis.
      virtual void SetDecoElement(int i, double val);
      /// Set mixing angle between two decoherence parameters a_i, a_j.
      virtual void SetDecoAngle(int i, int j, double th);

      /// Get the currently set power-law index for decoherence parameters.
      virtual int GetPower();
      /// Get the value of the a_i decoherence parameter in Gell-Mann basis.
      virtual double GetDecoElement(int i);
      /// Get mixing angle between two decoherence parameters a_i, a_j.
      virtual double GetDecoAngle(int i, int j);

      /// Get element i,j of the dissipator matrix in Gell-Mann basis.
      virtual double GetHGM(int i, int j);
      /// Get dissipator element in Gell-Mann basis.
      virtual double GetDissipatorElement(int i, int j);

      /// Reimplemented from PMNS_Base
      virtual void SetIsNuBar(bool isNuBar);
      /// Reimplemented from PMNS_DensityMatrix
      virtual matrixD ProbMatrix(int nflvi, int nflvf);

    protected:
      /// Reimplemented from PMNS_Base.
      virtual void BuildHms();
      /// Build effective Hamiltonian in Vacuum Mass Basis (VMB).
      virtual void BuildHVMB(NuPath p);

      /// Build effective Hamiltonian in Gell-Mann representation (GM).
      virtual void BuildHGM(NuPath p);
      /// Build the dissipator in Gell-Mann representation.
      virtual void BuildDissipator();
      /// Build the matrix M used in evolution equations.
      virtual void BuildM(NuPath p);

      /// Rotate rho to/from mass basis.
      virtual void RotateState(bool to_mass);
      /// Build Gell-Mann representation of density matrix.
      virtual void BuildR();
      /// Update density matrix from Gell-Mann representation.
      virtual void UpdateRho();

      /// Reimplemented from PMNS_Base
      virtual void Propagate();
      /// Reimplemented from PMNS_Base
      virtual void PropagatePath(NuPath p);

      int     fPower; ///< Power-law index (n).
      vectorD fa;     ///< Dissipator parameters |a_i|
      matrixD fcos;   ///< cosines ai . aj

      vectorD fR;    ///< Initial state of density matrix in Gell-Mann basis.
      vectorD fRt;   ///< Time evolution of density matrix in Gell-Mann basis.
      matrixC fUM;   ///< PMNS Matrix
      matrixC fHVMB; ///< Effective Hamiltonian in VMB (3x3).
      matrixD fHGM;  ///< Effective Hamiltonian in GMB (9x9).
      matrixD fD;    ///< Off-diagonal, 9x9 dissipator

      Eigen::MatrixXd fM; ///< Buffer matrix for exponential

      bool fBuiltDissipator; ///< Flag to rebuild dissipator
  };

} // namespace OscProb

#endif

///////////////////////////////////////////////////////////////////////////////
