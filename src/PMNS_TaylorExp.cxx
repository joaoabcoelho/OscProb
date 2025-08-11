///////////////////////////////////////////////////////////////////////////////
//
// Implementation of oscillations of neutrinos in matter in a
// three-neutrino framework with a first order Taylor expansion.
//
// This  class inherits from the PMNS_Fast class
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include "PMNS_TaylorExp.h"
#include "MatrixDecomp/zheevh3.h"
#include "PremModel.h"
#include <algorithm>


using namespace OscProb;

using namespace std;

//.............................................................................
///
/// Constructor. \sa PMNS_Base::PMNS_Base
///
/// This class is restricted to 3 neutrino flavours.
///
PMNS_TaylorExp::PMNS_TaylorExp() : PMNS_Fast() 
{
    InitializeTaylorsVectors();

    SetwidthBin(0,0);

    SetAvgProbPrec(1e-4);
}

//.............................................................................
///
/// Nothing to clean.
///
PMNS_TaylorExp::~PMNS_TaylorExp() {}

//.............................................................................
///
/// Set vector sizes and initialize elements to zero.
/// Initialize diagonal elements of S to one
///
void PMNS_TaylorExp::InitializeTaylorsVectors()
{
    matrixC densityMatrix = matrixC(fNumNus, vectorC(fNumNus, 0));

    flambdaInvE = vectorD(fNumNus, 0);
    fVInvE = matrixC(fNumNus, vectorC(fNumNus, 0));

    flambdaCosT = vectorD(fNumNus, 0);
    fVcosT = matrixC(fNumNus, vectorC(fNumNus, 0));

    fevolutionMatrixS = matrixC(fNumNus, vectorC(fNumNus, 0));

    fSflavor= matrixC(fNumNus, vectorC(fNumNus, 0));
    Kmass= matrixC(fNumNus, vectorC(fNumNus, 0));
    Kflavor= matrixC(fNumNus, vectorC(fNumNus, 0));

    for(int i= 0 ; i<fNumNus; i++){

        fevolutionMatrixS[i][i] = 1;

        for(int j = 0 ; j<fNumNus ; j++){
            fKInvE[i][j] = 0;
            fKcosT[i][j] = 0;
        }
    }
}

//.............................................................................
///
/// Set neutrino angle.
///
/// @param cosT - The cosine of the neutrino angle
///
void PMNS_TaylorExp::SetCosT(double cosT)
 {
    fcosT = cosT;
 }

//.............................................................................
///
/// Set bin's widths.
///
/// @param dE - The width of the bin for energy in GeV
/// @param dcosT - The width of the bin for angle
///
void PMNS_TaylorExp::SetwidthBin(double dE , double dcosT)
{
    fdInvE = dE;
    fdcosT = dcosT;
}


//.............................................................................
///
/// Build K matrix for the inverse of energy in mass basis. 
///
/// The variable for which a Taylor expansion is done here is not directly the 
/// energy but the inverse of it. This change of variable allow to express the
/// hamiltonien as linear with respect to this new variable.
///  
/// @param L - The length of the layer in GeV-1
/// @param K - The K matrix for the inverse of energy in the mass basis
///
void PMNS_TaylorExp::BuildKE(double L , matrixC& K)
{
    double lenghtEV = L * kKm2eV ; // L in eV-1
    double bufK =  lenghtEV * 0.5 ; // L/2 in eV-1

    for(int i = 0 ; i<fNumNus ; i++){

        complexD buffer[3];

        complexD Hms_kl;

        for(int l = 0 ; l<fNumNus ; l++){
            for(int k = 0 ; k<fNumNus ; k++){
                
                if (k<=l)
                    Hms_kl  = fHms[k][l];
                else
                    Hms_kl  = conj(fHms[l][k]);

                if(fIsNuBar && k!=l)
                    Hms_kl = conj(Hms_kl);

                buffer[l] += conj(fEvec[k][i]) *  Hms_kl;
            }
        }

        for(int j = 0 ; j<=i ; j++){

            K[i][j] = 0;

            for (int l = 0 ; l<fNumNus ; l++){
                K[i][j] += buffer[l] * fEvec[l][j] ;
            }
            
            complexD C;

            if(i == j){
                C = complexD(1,0);
            }
            else {
                double arg = (fEval[i] - fEval[j]) * lenghtEV ;

                C = complexD(1,0) * ( exp(complexD(0.0 , arg)) - complexD(1 , 0.0) ) / complexD(0.0 , arg)  ;
                
            }  

            K[i][j] *= bufK * C ;

            if(i != j)  
                K[j][i] = conj(K[i][j]);
            
        }
    }

}

//.............................................................................
///
/// Build K matrix for angle in flavor basis 
///  
/// @param L - The length of the layer in GeV-1
/// @param K - The K matrix for the angle in the flavor basis
///
void PMNS_TaylorExp::BuildKcosT(double L, matrixC& K)
{
    double theta = acos(fcosT);

    //cout<<fcosT<<"   "<<sin(theta)<<endl;


    for(int j = 0 ; j<fNumNus ; j++){
        for(int i = 0 ; i<=j ; i++){
            K[i][j] = ( -2* 6371 * kKm2eV * sin(theta) ) * fHam[i][j]; //  abs()???  
            //K[i][j] =   kKm2eV *  2 * L * fHam[i][j];

            if(i != j){
                K[j][i] = conj(K[i][j]);
            }
        }
    } 

    //printMatrix1(K);
}

//.............................................................................
///
/// Rotate the S matrix from mass to flavor basis
///
/// @param fPhases - The diagonal elements of S in mass basis
/// @param S - The S matrix in flavor basis
///
void PMNS_TaylorExp::rotateS(vectorC fPhases, matrixC& S)
{
    complexD buffer[3];

    for(int j = 0 ; j<fNumNus ; j++){ 

        for(int k = 0 ; k<fNumNus ; k++) { buffer[k] = fPhases[k] * conj(fEvec[j][k]); }
        
        for(int i = 0 ; i<fNumNus ; i++){ 
            S[i][j] = 0;
            for(int k = 0 ; k<fNumNus ; k++){
                S[i][j] += fEvec[i][k] * buffer[k];
            }
        }
    }
}

//.............................................................................
///
/// Rotate the K matrix from mass to flavor basis
///
/// @param Kmass - The K matrix in mass basis
/// @param Kflavor - The K matrix in flavor basis
///
void PMNS_TaylorExp::rotateK(matrixC Kmass , matrixC& Kflavor)
{
    for(int j = 0 ; j<fNumNus ; j++){

        complexD buffer[3];

        for(int k = 0 ; k<fNumNus ; k++){
            for(int l = 0 ; l<fNumNus ; l++){
                buffer[k] += Kmass[k][l] * conj(fEvec[j][l]);
            }
        }

        for(int i = 0 ; i<=j ; i++){

            Kflavor[i][j] = 0;

            for(int k = 0 ; k<fNumNus ; k++){
                Kflavor[i][j] += fEvec[i][k] * buffer[k];
            }

            if(i != j){ Kflavor[j][i] = conj(Kflavor[i][j]); }
        }
    }    
    
}


//.............................................................................
///
/// Product between two S matrices. 
///
/// This is used to calculate the matrix S corresponding to the propagation 
/// between the beginning of the path and the end of the current layer. 
///
/// The matrix fevolutionMatrixS represent the propagation between the beginning 
/// of the path and the beginning of the current layer. This matrix is updated 
/// after every layer with this function.
///
/// @param SLayer - The S matrix corresponding to the propagation in the current 
///                 layer
///
void PMNS_TaylorExp::MultiplicationRuleS(matrixC SLayer)
{
    complexD save [3];

    for(int j = 0 ; j<fNumNus ; j++){

        for(int n = 0 ; n <fNumNus ; n++) { save[n] = fevolutionMatrixS[n][j]; }

        for(int i = 0 ; i<fNumNus ; i++){ 

            fevolutionMatrixS[i][j] = 0;

            for(int k = 0 ; k<fNumNus ; k++){
                fevolutionMatrixS[i][j] += SLayer[i][k] * save[k];
            }
        }
    }

}

//.............................................................................
///
/// Product between two K matrices. 
///
/// This is used to calculate the matrix K corresponding to the propagation 
/// between the beginning of the path and the end of the current layer. 
///
/// @param KLayer - The S matrix corresponding to the propagation in the current 
///                 layer
/// @param K - The S matrix corresponding to the propagation between the beginning 
///            of the path and the beginning of the current layer
///
void PMNS_TaylorExp::MultiplicationRuleK(matrixC KLayer,complexD K[3][3])
{

    for(int i = 0 ; i<fNumNus ; i++){

        complexD buffer[3];

        for(int l = 0 ; l<fNumNus ; l++){
            for(int k = 0 ; k<fNumNus ; k++){
                buffer[l] += conj(fevolutionMatrixS[k][i]) * KLayer[k][l];
            }
        }

        for(int j = 0 ; j<=i ; j++){

            for(int l = 0 ; l<fNumNus ; l++){                                                     
                K[i][j] += buffer[l] * fevolutionMatrixS[l][j] ; 
            }
            
            if(i != j){
                K[j][i] = conj(K[i][j]); 
            }
        }
    }  

}


//.............................................................................
///
/// Propagate neutrino state through full path
///
void PMNS_TaylorExp::PropagateTaylor()
{
  for (int i = 0; i < int(fNuPaths.size()); i++) { PropagatePathTaylor(fNuPaths[i]); }
}

//.............................................................................
///
/// Propagate the current neutrino state through a given path
/// @param p - A neutrino path segment
///
void PMNS_TaylorExp::PropagatePathTaylor(NuPath p)
{
    // Set the neutrino path
    SetCurPath(p);

    // Solve for eigensystem
    SolveHam();   

    // Get the evolution matrix in mass basis
    double LengthIneV = kKm2eV * p.length;      
    for (int i = 0; i < fNumNus; i++) {
        double arg = fEval[i] * LengthIneV;
        fPhases[i] = complexD(cos(arg), -sin(arg));
    }

    // Rotate S in flavor basis
    //matrixC fSflavor = matrixC(fNumNus, vectorC(fNumNus, 0));
    rotateS(fPhases,fSflavor);                                    

    // if avg on E
    if(fdInvE != 0){

        // Build KE in mass basis
        //matrixC Kmass = matrixC(fNumNus, vectorC(fNumNus, 0));
        BuildKE(p.length,Kmass);

        // Rotate KE in flavor basis
        //matrixC Kflavor = matrixC(fNumNus, vectorC(fNumNus, 0));
        rotateK(Kmass,Kflavor);

        // Multiply this layer K's with the previous path K's
        MultiplicationRuleK(Kflavor,fKInvE);    
        
    }

    // if avg on cosT
    if(fdcosT != 0){
        matrixC Kmass = matrixC(fNumNus, vectorC(fNumNus, 0));
        BuildKcosT(p.length, Kmass);

        // Multiply this layer K's with the previous path K's
        MultiplicationRuleK(Kmass,fKcosT);

        //********************************************************************************************* */

        //fCountLayer++;

        //*********************************************************************************************
        
    }

    // Multiply this layer S's with the previous path S's
    MultiplicationRuleS(fSflavor);

}

//.............................................................................
///
/// Solve one K matrix for eigenvectors and eigenvalues.
///
/// This is using a method from the GLoBES software available at
/// http://www.mpi-hd.mpg.de/personalhomes/globes/3x3/ \n
/// We should cite them accordingly
///
/// @param K - The K matrix 
/// @param lambda - The eigenvalues of K 
/// @param V - The eigenvectors of K 
///
void PMNS_TaylorExp::SolveK(complexD K[3][3], vectorD& lambda, matrixC& V)
{
    double   fEvalGLoBES[3];
    complexD fEvecGLoBES[3][3];

    // Solve K for eigensystem using the GLoBES method
    zheevh3(K, fEvecGLoBES, fEvalGLoBES);

    // Fill flambdaInvE and fVInvE vectors from GLoBES arrays
    for (int i = 0; i < fNumNus; i++) {
        lambda[i] = fEvalGLoBES[i]; 
        for (int j = 0; j < fNumNus; j++) { V[i][j] = fEvecGLoBES[i][j]; }
    }
}

//.............................................................................
///
/// Formula for the average probability of flvi going to flvf over 
/// a bin of energy E with width dE with a Taylor expansion.
///
/// Flavours are:
/// <pre>
///   0 = nue, 1 = numu, 2 = nutau
///   3 = sterile_1, 4 = sterile_2, etc.
/// </pre>
///
/// @param flvi - The neutrino starting flavour.
/// @param flvf - The neutrino final flavour.
/// @param dbin - The width of the bin
/// @param lambda - The eigenvalues of K 
/// @param V - The eigenvectors of K 
///
/// @return Average neutrino oscillation probability
///
double PMNS_TaylorExp::AvgFormula(int flvi, int flvf, double dbin, vectorD lambda, matrixC V)
{

    vectorC SV = vectorC(fNumNus, 0);

    for(int i = 0 ; i<fNumNus ; i++) {
        for(int k = 0 ; k<fNumNus ; k++){
            SV[i] += fevolutionMatrixS[flvf][k] *V[k][i]; 
        }
    }
    
    complexD buffer [3];

    for( int n = 0 ; n<fNumNus ; n++){
        buffer[n] = SV[n] * conj(V[flvi][n]);
    }

    complexD sinc[fNumNus][fNumNus];

    for(int j = 0 ; j<fNumNus ; j++){

        sinc[j][j] = 1;
        
        for(int i = 0 ; i<j ; i++){
            double arg = (lambda[i] - lambda[j]) * dbin / 2; 
            sinc[i][j] = sin(arg) / arg;
            sinc[j][i] = sinc[i][j];
        }
    }
    
    complexD P = 0;

    for(int j = 0 ; j<fNumNus ; j++){
        for(int i = 0 ;i<fNumNus ;i++){
            P += buffer[i] *conj(buffer[j]) * sinc[j][i];
        }
    }

    return real(P); 
}

//.............................................................................
///
/// Compute the average probability of flvi going to flvf over 
/// a bin of energy E with width dE with a Taylor expansion.
///
/// This gets transformed into L/E, since the oscillation terms
/// have arguments linear in 1/E and not E. 
///
/// Flavours are:
/// <pre>
///   0 = nue, 1 = numu, 2 = nutau
///   3 = sterile_1, 4 = sterile_2, etc.
/// </pre>
///
/// @param flvi - The neutrino starting flavour.
/// @param flvf - The neutrino final flavour.
/// @param E - The neutrino energy in GeV
/// @param dE - The energy bin width in GeV
///
/// @return Average neutrino oscillation probability
///
double PMNS_TaylorExp::AvgProb(int flvi, int flvf, double E , double dE)
{
    if (E <= 0) return 0;

    if (fNuPaths.empty()) return 0;

    // Don't average zero width
    if (dE <= 0) return Prob(flvi, flvf, E);

    vectorD Ebin = ConvertEtoLoE(E,dE);

    //cout<<endl<<"(NEW BIN)-------------------------E = "<<E<<"-----------------------------------"<<endl;

    //return fct avr proba
    return AvgProbLoE(flvi, flvf, Ebin[0], Ebin[1]);
}

//.............................................................................
///
/// Compute the average probability of flvi going to flvf over 
/// a bin of energy L/E with width dLoE with a Taylor expansion.
///
/// Flavours are:
/// <pre>
///   0 = nue, 1 = numu, 2 = nutau
///   3 = sterile_1, 4 = sterile_2, etc.
/// </pre>
///
/// @param flvi - The neutrino starting flavour.
/// @param flvf - The neutrino final flavour.
/// @param LoE - The neutrino  L/E value in the bin center in km/GeV
/// @param dLoE - The L/E bin width in km/GeV
///
/// @return Average neutrino oscillation probability
///
double PMNS_TaylorExp::AvgProbLoE(int flvi, int flvf, double LoE , double dLoE)
{
    if (LoE <= 0) return 0;

    if (fNuPaths.empty()) return 0;

    // Don't average zero width
    if (dLoE <= 0) return Prob(flvi, flvf, fPath.length / LoE);

    // Get sample points for this bin
    vectorD samples = GetSamplePoints(LoE, dLoE);

    double avgprob  = 0;
    double L = fPath.length;
    double sumw   = 0;

    // Loop over all sample points
    for (int j = 1; j < int(samples.size()); j++) {

        double w = 1. / pow(samples[j], 2);

        avgprob += w * AvgAlgo(flvi, flvf, samples[j], samples[0],L);

        sumw += w;

    }
    
    // Return average of probabilities
    return  avgprob / sumw;
}

//.............................................................................
///
/// 
///
double PMNS_TaylorExp::AvgAlgo(int flvi, int flvf, double LoE , double dLoE, double L)
{
    // Set width bin as 1/E 
    //double L = fPath.length;
    double d1oE = dLoE / L;

    // reset K et S et Ve et lambdaE
    InitializeTaylorsVectors();

    SetEnergy(L / LoE);
    SetwidthBin(d1oE, 0);

    //Propagate -> get S and K matrix (on the whole path)
    PropagateTaylor();

    //DiagolK -> get VE and lambdaE
    SolveK(fKInvE,flambdaInvE,fVInvE);

    //return fct avr proba
    return AvgFormula(flvi, flvf, d1oE/ kGeV2eV , flambdaInvE, fVInvE);
}

//.............................................................................
///
/// Compute the average probability of flvi going to flvf over 
/// a bin of angle cost with width dcosT with a Taylor expansion.
///
/// Flavours are:
/// <pre>
///   0 = nue, 1 = numu, 2 = nutau
///   3 = sterile_1, 4 = sterile_2, etc.
/// </pre>
///
/// @param flvi - The neutrino starting flavour.
/// @param flvf - The neutrino final flavour.
/// @param E - The neutrino energy in GeV
/// @param cosT - The cosine of the neutrino angle
/// @param dcosT - The cosT bin width in GeV-1
///
/// @return Average neutrino oscillation probability
///
double PMNS_TaylorExp::AvgProbAngle(int flvi, int flvf, double E, double cosT , double dcosT)
{
    // reset K et S et Ve et lambdaE
    InitializeTaylorsVectors();

    SetEnergy(E);
    SetCosT(cosT);
    SetwidthBin(0,dcosT);

    //Propagate -> get S and K matrix (on the whole path)
    PropagateTaylor();

    //DiagolK -> get VE and lambdaE
    SolveK(fKcosT,flambdaCosT,fVcosT);

    //return fct avr proba
    return AvgFormula(flvi,flvf,fdcosT, flambdaCosT, fVcosT);
}

//.............................................................................
///
/// Fomula for the propability for flvi going to flvf for an energy E+dE 
/// using a first order Taylor expansion from a reference energy E.
///
/// Flavours are:
/// <pre>
///   0 = nue, 1 = numu, 2 = nutau
///   3 = sterile_1, 4 = sterile_2, etc.
/// </pre>
///
/// @param flvi - The neutrino starting flavour.
/// @param flvf - The neutrino final flavour.
/// @param dE - The energy variation in GeV
/// @param lambda - The eigenvalues of K 
/// @param V - The eigenvectors of K 
///
/// @return Neutrino oscillation probability
///
double PMNS_TaylorExp::AvgFormulaExtrapolation(int flvi, int flvf, double dbin, vectorD lambda, matrixC V)
{
 
    vectorC SV = vectorC(fNumNus, 0);

    for(int i = 0 ; i<fNumNus ; i++) {
        for(int k = 0 ; k<fNumNus ; k++){
            SV[i] += fevolutionMatrixS[flvf][k] *V[k][i]; 
        }
    }
    
    complexD buffer [3];

    for( int n = 0 ; n<fNumNus ; n++){
        buffer[n] = SV[n] * conj(V[flvi][n]);
    }

    complexD expo[fNumNus][fNumNus];

    for(int j = 0 ; j<fNumNus ; j++){
        
        for(int i = 0 ; i<fNumNus ; i++){
            double arg = (lambda[j] - lambda[i]) * dbin ;    
            expo[j][i] = exp(complexD(0.0, arg)) ;
        }
    }
    
    complexD P = 0;

    for(int j = 0 ; j<fNumNus ; j++){
        for(int i = 0 ;i<fNumNus ;i++){
            P += buffer[i] *conj(buffer[j]) * expo[j][i];
        }
    }

    return real(P);

}

//.............................................................................
///
/// Compute the propability for flvi going to flvf for an energy E+dE 
/// using a first order Taylor expansion from a reference energy E.
///
/// Flavours are:
/// <pre>
///   0 = nue, 1 = numu, 2 = nutau
///   3 = sterile_1, 4 = sterile_2, etc.
/// </pre>
///
/// @param flvi - The neutrino starting flavour.
/// @param flvf - The neutrino final flavour.
/// @param E - The reference energy in GeV
/// @param dE - The energy variation in GeV
///
/// @return Neutrino oscillation probability
///
double PMNS_TaylorExp::interpolationEnergy(int flvi, int flvf, double E , double dE)
{
    // reset K et S et Ve et lambdaE
    InitializeTaylorsVectors();

    SetEnergy(E);
    SetwidthBin( 1 / (E + dE) - 1 / E ,0);

    //Propagate -> get S and K matrix (on the whole path)
    PropagateTaylor();

    //DiagolK -> get VE and lambdaE
    SolveK(fKInvE,flambdaInvE,fVInvE);

    return AvgFormulaExtrapolation(flvi,flvf,fdInvE/kGeV2eV, flambdaInvE, fVInvE);
}

//.............................................................................
///
///
///
double PMNS_TaylorExp::interpolationEnergyLoE(int flvi, int flvf, double LoE , double dLoE)
{
    // reset K et S et Ve et lambdaE
    InitializeTaylorsVectors();

    SetCurPath(AvgPath(fNuPaths));
    double L = fPath.length;

    SetEnergy(LoE);
    SetwidthBin(dLoE,0);

    //Propagate -> get S and K matrix (on the whole path)
    PropagateTaylor();

    //DiagolK -> get VE and lambdaE
    SolveK(fKInvE,flambdaInvE,fVInvE);

    return AvgFormulaExtrapolation(flvi,flvf,dLoE*kGeV2eV, flambdaInvE, fVInvE);
}


//.............................................................................
///
/// Compute the propability for flvi going to flvf for an angle cosT+dcosT 
/// using a first order Taylor expansion from a reference angle cosT.
///
/// Flavours are:
/// <pre>
///   0 = nue, 1 = numu, 2 = nutau
///   3 = sterile_1, 4 = sterile_2, etc.
/// </pre>
///
/// @param flvi - The neutrino starting flavour.
/// @param flvf - The neutrino final flavour.
/// @param cosT - The reference angle 
/// @param dcosT - The angle variation 
///
/// @return Neutrino oscillation probability
///
double PMNS_TaylorExp::interpolationCosT(int flvi, int flvf, double cosT , double dcosT)
{
    // reset K et S et Ve et lambdaE
    InitializeTaylorsVectors();

    //SetEnergy(E);
    SetCosT(cosT);
    SetwidthBin(0,dcosT);

    //Propagate -> get S and K matrix (on the whole path)
    PropagateTaylor();

    //DiagolK -> get VE and lambdaE
    SolveK(fKcosT,flambdaCosT,fVcosT);

    return AvgFormulaExtrapolation(flvi,flvf,fdcosT, flambdaCosT, fVcosT);
}





//.............................................................................
///
/// Compute the sample points for a bin of L/E with width dLoE
///
/// This is used for averaging the probability over a bin of L/E..
///
/// @param LoE  - The neutrino  L/E value in the bin center in km/GeV
/// @param dLoE   - The L/E bin width in km/GeV
///
vectorD PMNS_TaylorExp::GetSamplePoints(double LoE, double dLoE)
{

  // Set a number of sub-divisions to achieve "good" accuracy
  // This needs to be studied better
  int n_div = ceil(10* pow(dLoE / LoE, 0.8) * pow(LoE, 0.2) /  sqrt(fAvgProbPrec / 1e-4));
  //int n_div = ceil(10* pow(dLoE / LoE, 0.9) * pow(LoE, 0.1) * (1 + 2 * exp(-pow((LoE-2100),2)/(2*pow(300,2))) ) /  sqrt(fAvgProbPrec / 1e-4));
  //int n_div = 15;
  //cout<<"nbr sub-bins = "<<n_div<<endl;

  // A vector to store sample points
  vectorD allSamples;
  allSamples.push_back(dLoE / n_div);

  // Loop over sub-divisions
  for (int k = 0; k < n_div; k++) {
    // Define sub-division center and width
    double bctr = LoE - dLoE / 2 + (k + 0.5) * dLoE / n_div;
   
    allSamples.push_back(bctr);

  } // End of loop over sub-divisions

  // Return all sample points
  return allSamples;
}








/*double PMNS_TaylorExp::AvgProb(int flvi, int flvf, double E , double dE, double a, double b, double c, double k)
{
    if (E <= 0) return 0;

    if (fNuPaths.empty()) return 0;

    // Don't average zero width
    if (dE <= 0) return Prob(flvi, flvf, E);

    vectorD Ebin = ConvertEtoLoE(E,dE);

    //cout<<endl<<"(NEW BIN)-------------------------E = "<<E<<"-----------------------------------"<<endl;

    //return fct avr proba
    return AvgProbLoE(flvi, flvf, Ebin[0], Ebin[1],a,b,c,k);
}

double PMNS_TaylorExp::AvgProbLoE(int flvi, int flvf, double LoE , double dLoE, double a, double b, double c, double k)
{
    if (LoE <= 0) return 0;

    if (fNuPaths.empty()) return 0;

    // Don't average zero width
    if (dLoE <= 0) return Prob(flvi, flvf, fPath.length / LoE);

    // Get sample points for this bin
    vectorD samples = GetSamplePoints(LoE, dLoE,a,b,c,k);
    
    double avgprob  = 0;
    double L = fPath.length;
    double sumw   = 0;

    // Loop over all sample points
    for (int j = 1; j < int(samples.size()); j++) {

        double w = 1. / pow(samples[j], 2);

        avgprob += w * AvgAlgo(flvi, flvf, samples[j], samples[0],L);

        sumw += w;

    }
    
    // Return average of probabilities
    return  avgprob / sumw;
}



vectorD PMNS_TaylorExp::GetSamplePoints(double LoE, double dLoE, double a, double b, double c, double k)
{

  // Set a number of sub-divisions to achieve "good" accuracy
  // This needs to be studied better
  int n_div = ceil(k* pow(dLoE / LoE, a) * pow(LoE, b) * (1 + 2 * exp(-pow((LoE-2100),2)/(2*pow(500,2))) ) / pow(fAvgProbPrec / 1e-4 , c));
  //int n_div = 15;
  //cout<<"nbr sub-bins = "<<n_div<<endl;

  // A vector to store sample points
  vectorD allSamples;
  allSamples.push_back(dLoE / n_div);

  // Loop over sub-divisions
  for (int k = 0; k < n_div; k++) {
    // Define sub-division center and width
    double bctr = LoE - dLoE / 2 + (k + 0.5) * dLoE / n_div;
   
    allSamples.push_back(bctr);

  } // End of loop over sub-divisions

  // Return all sample points
  return allSamples;
}*/




//.............................................................................
///
///
///
void PMNS_TaylorExp::RotateDensityM(bool to_mass, matrixC V, matrixC& densityMatrix)
{
    matrixC Buffer = matrixC(fNumNus, vectorC(fNumNus, 0));

  // buffer = rho . U
  for (int i = 0; i < fNumNus; i++) {
    for (int j = 0; j < fNumNus; j++) {
      for (int k = 0; k < fNumNus; k++) {
        if (to_mass)
          Buffer[i][j] += densityMatrix[i][k] * V[k][j];
        else
          Buffer[i][j] += densityMatrix[i][k] * conj(V[j][k]);
      }
    }
  }

  // rho = U^\dagger . buffer = U^\dagger . rho . U
  // Final matrix is Hermitian, so copy upper to lower triangle
  for (int i = 0; i < fNumNus; i++) {
    for (int j = i; j < fNumNus; j++) {
        densityMatrix[i][j] = 0;
      for (int k = 0; k < fNumNus; k++) {
        if (to_mass)
            densityMatrix[i][j] += conj(V[k][i]) * Buffer[k][j];
        else
            densityMatrix[i][j] += V[i][k] * Buffer[k][j];
      }
      if (j > i) densityMatrix[j][i] = conj(densityMatrix[i][j]);
    }
  }
}

//.............................................................................
///
///
///
void PMNS_TaylorExp::HadamardProduct(vectorD lambda, matrixC& densityMatrix, double dbin)
{
    matrixC sinc = matrixC(fNumNus, vectorC(fNumNus, 0));
    for(int j=0 ; j<fNumNus ; j++){
        for(int i = 0 ; i<j ; i++){
            double arg = (lambda[i] - lambda[j]) * dbin * kGeV2eV; //ATTTTTTTTTTTTTTTTTTTTTTT
            sinc[i][j] = sin(arg)/arg;

            sinc[j][i] = sinc[i][j];
        }
    }    // PEUT ETRE AMELIORE !!!!!!!!!!!!!!!!!!!!!!

    for(int i=0 ; i<fNumNus ; i++){
        sinc[i][i] = 1;
    }
    
    for(int j=0 ; j<fNumNus ; j++){
        for(int i = 0 ; i<fNumNus ; i++){
            densityMatrix[i][j] = densityMatrix[i][j] * sinc[i][j];
        }
    }
}

//.............................................................................
///
///
///
double PMNS_TaylorExp::AvgAlgorithm(int flvi, int flvf)
{
    //matrixC densityMatrix = matrixC(fNumNus, vectorC(fNumNus, 0));
    densityMatrix[flvi][flvi] = 1;

    RotateDensityM(true,fVcosT,densityMatrix);
    HadamardProduct(flambdaCosT,densityMatrix,fdcosT);
    RotateDensityM(false,fVcosT,densityMatrix);

    RotateDensityM(true,fVcosT,densityMatrix);
    HadamardProduct(flambdaInvE,densityMatrix,fdInvE);
    RotateDensityM(false,fVcosT,densityMatrix);

    RotateDensityM(false,fevolutionMatrixS,densityMatrix);

    return real(densityMatrix[flvf][flvf]);
}

//.............................................................................
///
///
///
double PMNS_TaylorExp::AvgProb(int flvi, int flvf, double E , double dE, double cosT , double dcosT)
{
    // reset K et S et Ve et lambdaE
    InitializeTaylorsVectors();

    SetEnergy(E);
    SetCosT(cosT);
    SetwidthBin(dE,dcosT);

    //Propagate -> get S and K matrix (on the whole path)
    PropagateTaylor();

    //DiagolK -> get VE and lambdaE
    SolveK(fKInvE,flambdaInvE,fVInvE);
    SolveK(fKcosT,flambdaCosT,fVcosT);

    return AvgAlgorithm(flvi,flvf);
}






















void PMNS_TaylorExp::printMatrix1(matrixC M)
{
    cout<<endl;
    for(int j=0 ; j<3 ; j++)
    {
        for(int k=0 ; k<3 ; k++)
        {
            cout<<M[j][k]<<" ";
        }
        cout<<endl;
    }
    cout<<endl;
}
void PMNS_TaylorExp::printMatrix2(complexD M[3][3])
{
    cout<<endl;
    for(int j=0 ; j<3 ; j++)
    {
        for(int k=0 ; k<3 ; k++)
        {
            cout<<M[j][k]<<" ";
        }
        cout<<endl;
    }
    cout<<endl;
}
void PMNS_TaylorExp::LenghtLayer()
{
    for(int i = 0 ; i<fNuPaths.size() ; i++)
    {
        cout<<"Layer "<<i<<"    L = "<<fNuPaths[i].length<<"   Density = "<<fNuPaths[i].density<<endl;
    }
}