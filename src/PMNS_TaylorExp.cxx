///////////////////////////////////////////////////////////////////////////////
//info
///////////////////////////////////////////////////////////////////////////////

#include <iostream>

#include "PMNS_TaylorExp.h"

#include "MatrixDecomp/zheevh3.h"

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
}

//.............................................................................
///
/// Nothing to clean.
///
PMNS_TaylorExp::~PMNS_TaylorExp() {}

//.............................................................................
///
///
///
void PMNS_TaylorExp::printMatrix1(matrixC M)
{
    for(int j=0 ; j<3 ; j++)
    {
        for(int k=0 ; k<3 ; k++)
        {
            cout<<M[j][k]<<" ";
        }
        cout<<endl;
    }
}
void PMNS_TaylorExp::printMatrix2(complexD M[3][3])
{
    for(int j=0 ; j<3 ; j++)
    {
        for(int k=0 ; k<3 ; k++)
        {
            cout<<M[j][k]<<" ";
        }
        cout<<endl;
    }
}

//.............................................................................
///
///
///
void PMNS_TaylorExp::InitializeTaylorsVectors()
{
    flambdaE = vectorD(fNumNus, 0);
    fVE = matrixC(fNumNus, vectorC(fNumNus, 0));

    flambdaCosT = vectorD(fNumNus, 0);
    fVcosT = matrixC(fNumNus, vectorC(fNumNus, 0));

    fevolutionMatrixS = matrixC(fNumNus, vectorC(fNumNus, 0));

    for(int i= 0 ; i<fNumNus; i++){

        fevolutionMatrixS[i][i] = 1;

        for(int j = 0 ; j<fNumNus ; j++){
            fKE[i][j] = 0;
            fKcosT[i][j] = 0;
        }
    }

}

void PMNS_TaylorExp::SetCosT(double cosT)
 {
    fcosT = cosT;
 }

//.............................................................................
///
///
///
void PMNS_TaylorExp::SetwidthBin(double dE , double dcosT)
{
    fdE = dE;
    fdcosT = dcosT;
}

//.............................................................................
///
///
///
void PMNS_TaylorExp::rotateS(vectorC fPhases,matrixC& S)
{
    for(int j = 0 ; j<fNumNus ; j++){ 

        complexD buffer[3];

        for(int k = 0 ; k<fNumNus ; k++){
            buffer[k] = fPhases[k] * conj(fEvec[j][k]);
        }
        
        for(int i = 0 ; i<fNumNus ; i++){ 
            for(int k = 0 ; k<fNumNus ; k++){
                S[i][j] += fEvec[i][k] * buffer[k];
            }
        }
    }
    /*
    for(int j = 0 ; j<fNumNus ; j++){ //PEUT ETRE AMELIORER
        for(int i = 0 ; i<fNumNus ; i++){ //CHG
            for(int k = 0 ; k<fNumNus ; k++){
                S[i][j] += fEvec[i][k] * fPhases[k] * conj(fEvec[j][k]);
            }
        }
    }*/
}

//.............................................................................
///
///
///
void PMNS_TaylorExp::rotateK(matrixC Kmass,matrixC& Kflavor)
{
    for(int j = 0 ; j<fNumNus ; j++){

        complexD par[3];

        for(int k = 0 ; k<fNumNus ; k++){
            for(int l = 0 ; l<fNumNus ; l++){
                par[k] += Kmass[k][l] * conj(fEvec[j][l]);
            }
        }

        for(int i = 0 ; i<=j ; i++){
            for(int k = 0 ; k<fNumNus ; k++){
                Kflavor[i][j] += fEvec[i][k] * par[k];
            }

            if(i != j){
                Kflavor[j][i] = conj(Kflavor[i][j]);
            }
        }
    }
    
    /*for(int j = 0 ; j<fNumNus ; j++)
    {
        for(int i = 0 ; i<=j ; i++)
        {
            for(int k = 0 ; k<fNumNus ; k++)
            {
                for(int l = 0 ; l<fNumNus ; l++)
                {
                    Kflavor[i][j] += fEvec[i][k] * Kmass[k][l] * conj(fEvec[j][l]);
                }
            }

            if(i != j){
                Kflavor[j][i] = conj(Kflavor[i][j]);
            }
        }
    }*/
    
    
}

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
void PMNS_TaylorExp::BuildKE(double L , matrixC& K)
{

    double lv =  kGeV2eV * fEnergy; // E in eV 

    for(int j = 0 ; j<fNumNus ; j++){
        for(int i = 0 ; i<=j ; i++){
            //K[i][j] = - kKm2eV * (L / (2*lv*lv)) * fHms[i][j]; 
            K[i][j] = - kKm2eV * (L / lv) * fHam[i][j];  


            if(i != j){
                K[j][i] = conj(K[i][j]);
            }
        }
    }

    /*double lv =  kGeV2eV * fEnergy; // E in eV 

    for(int j = 0 ; j<fNumNus ; j++){
        for(int i = 0 ; i<=j ; i++){
            for(int k = 0 ; k<fNumNus ; k++){
                for(int l = 0 ; l<fNumNus ; l++){
                    K[i][j] += conj(fEvec[k][i]) * fHms[k][l] * fEvec[l][j];
                }
            }
            K[i][j] *= - (1 / (2*lv*lv));

            complexD C;
            if(i == j){
                C = complexD(1,0);
            }
            else{
                double argg = (fEval[i] - fEval[j]) * L * kKm2eV;
                C = (complexD(cos(argg), sin(argg)) - complexD(1,0) ) / (complexD(0,argg)); 
            }// C=(1,0) because of H H' commutation (due to cst density case)

            K[i][j] *= kKm2eV * L  * K[i][j] * C; 

            if(i != j){
                K[j][i] = conj(K[i][j]);
            }

        }
    }*/

    /*matrixC Hprime = matrixC(fNumNus, vectorC(fNumNus, 0));
    for (int j = 0; j < fNumNus; j++) {
        // Set mass splitting
        Hprime[j][j] = fDm[j];
        // Reset off-diagonal elements
        for (int i = 0; i < j; i++) { Hprime[i][j] = 0; }
        // Rotate j neutrinos
        for (int i = 0; i < j; i++) { RotateH(i, j, Hprime); }
    }

    double lv =  kGeV2eV * fEnergy;

    for(int j = 0 ; j<fNumNus ; j++){
        for(int i = 0 ; i<=j ; i++){
            K[i][j] = - kKm2eV * (L / (lv*lv)) * Hprime[i][j];

            if(i != j){
                K[j][i] = conj(K[i][j]);
            }
        }
    }*/

}

//.............................................................................
///
///
///
void PMNS_TaylorExp::BuildKcosT(matrixC& K)
{
    double theta = acos(fcosT);

    for(int j = 0 ; j<fNumNus ; j++){
        for(int i = 0 ; i<=j ; i++){
            K[i][j] = 6371 * abs(sin(theta)) * fHam[i][j]; // 6371 en km ici

            if(i != j){
                K[j][i] = conj(K[i][j]);
            }
        }
    }
}

//.............................................................................
///
///
///
void PMNS_TaylorExp::MultiplicationRuleK(matrixC KLayer,complexD K[3][3])
{
    for(int j = 0 ; j<fNumNus ; j++){
        for(int i = 0 ; i<=j ; i++){

            for(int k = 0 ; k<fNumNus ; k++){
                for(int l = 0 ; l<fNumNus ; l++){                                                          
                    K[i][j] += conj(fevolutionMatrixS[k][i]) * KLayer[k][l] * fevolutionMatrixS[l][j] ; 
                }
            }
            
            if(i != j){
                K[j][i] = conj(K[i][j]); 
            }
        }
    }    

}

//.............................................................................
///
///
///
void PMNS_TaylorExp::MultiplicationRuleS(matrixC SLayer)
{
    matrixC SCopy = matrixC(fNumNus, vectorC(fNumNus, 0));
    for(int j = 0 ; j<fNumNus ; j++){
        for(int i = 0 ; i<fNumNus ; i++){
            SCopy[i][j] = fevolutionMatrixS[i][j];
        }
    }

    for(int j = 0 ; j<fNumNus ; j++){
        for(int i = 0 ; i<fNumNus ; i++){ 
            fevolutionMatrixS[i][j] = 0;
            for(int k = 0 ; k<fNumNus ; k++){

                fevolutionMatrixS[i][j] += SLayer[i][k] * SCopy[k][j];
            }
        }
    }

    //S
    /*
    for(int j = 0 ; j<fNumNus ; j++)
    {
        for(int i = 0 ; i<=j ; i++)
        {
            for(int k = 0 ; k<fNumNus ; k++)
            {
                if(j>k){
                    fevolutionMatrixS[i][j] -= SLayer[i][k] * conj(fevolutionMatrixS[j][k]);
                }
                else{
                    fevolutionMatrixS[i][j] += SLayer[i][k] * fevolutionMatrixS[k][j];
                }
            }
        }
    }
    //PEUT ETRE AMELIORER
    for(int j = 1 ; j<fNumNus ; j++)
    {
        for(int i = 0 ; i<j ; i++)
        {
            fevolutionMatrixS[j][i] = -conj(fevolutionMatrixS[i][j]);
        }
    }
    */
}

//.............................................................................
///
///
///
void PMNS_TaylorExp::LenghtLayer()
{
    for(int i = 0 ; i<fNuPaths.size() ; i++)
    {
        cout<<"Layer "<<i<<"    L = "<<fNuPaths[i].length<<"   Density = "<<fNuPaths[i].density<<endl;
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
    double LengthIneV = kKm2eV * p.length;      //IL FAUT CIONVERTIR TOUS LES L
    for (int i = 0; i < fNumNus; i++) {
        double arg = fEval[i] * LengthIneV;
        fPhases[i] = complexD(cos(arg), -sin(arg));
    }

    // Rotate S in flavor basis
    matrixC Sflavor = matrixC(fNumNus, vectorC(fNumNus, 0));
    rotateS(fPhases,Sflavor);                                    

    // Build KE in mass basis
    if(fdE != 0){
        matrixC Kmass = matrixC(fNumNus, vectorC(fNumNus, 0));
        BuildKE(p.length,Kmass);

        // Rotate KE in flavor basis
        //matrixC Kflavor = matrixC(fNumNus, vectorC(fNumNus, 0));
        //rotateK(Kmass,Kflavor);

        //multiplication rule for K and S 
        MultiplicationRuleK(Kmass,fKE);     //MARCHE POUR SYM CAR COMMUTTE AVEC H @@@@@@@@@@@@@@@@@@@@@@//
        
    }

    if(fdcosT != 0){
        matrixC Kmass = matrixC(fNumNus, vectorC(fNumNus, 0));
        BuildKcosT(Kmass);

        // Rotate KE in flavor basis
        //matrixC Kflavor = matrixC(fNumNus, vectorC(fNumNus, 0));
        //rotateK(Kmass,Kflavor);

        //multiplication rule for K and S 
        MultiplicationRuleK(Kmass,fKcosT);
        
    }

    MultiplicationRuleS(Sflavor);

}

//.............................................................................
///
/// 
///
void PMNS_TaylorExp::SolveK(complexD K[3][3], vectorD& lambda, matrixC& V)
{
    double   fEvalGLoBES[3];
    complexD fEvecGLoBES[3][3];

    // Solve Hamiltonian for eigensystem using the GLoBES method
    zheevh3(K, fEvecGLoBES, fEvalGLoBES);

    // Fill fEval and fEvec vectors from GLoBES arrays
    for (int i = 0; i < fNumNus; i++) {
        lambda[i] = fEvalGLoBES[i];
        for (int j = 0; j < fNumNus; j++) { V[i][j] = fEvecGLoBES[i][j]; }
    }

    //IL MANQUE DES CHOSES ICI
}

//.............................................................................
///
/// 
///
double PMNS_TaylorExp::avgFormula(int flvi, int flvf, double dbin, vectorD lambda, matrixC V)
{
    /*vectorC SVmulti = vectorC(fNumNus, 0);
    for(int i = 0 ; i<fNumNus ; i++) {
        for(int j = 0 ; j<fNumNus ; j++){
            SVmulti[i] += fevolutionMatrixS[flvf][j] *V[j][i]; //SOMMER SUR K?????????
        }
    }

    complexD P;
    complexD s1[fNumNus]; //le plus rapide? vector+pushback OU tableau+size definie des le débuts?
    //complexD s2[fNumNus];
    complexD sinc[fNumNus][fNumNus];

    for(int i = 0 ; i<fNumNus ; i++){
        s1[i] = SVmulti[i] * conj(V[flvi][i]) ;
        //s2[i] = conj(SVmulti[i]) * fV[flvi][i] ;  S1+CONJ(S2)!!!!!!!!!!!!!

        //sinc[i][i] = 1;
        for(int j = 0 ; j<i ; j++){
            double arg = (lambda[j] - lambda[i]) * dbin;
            sinc[j][i] = sin(arg)/arg;
            cout<<arg<<"    ";
        }
    }

    for(int j = 0 ; j<fNumNus ; j++){
        for(int i = 0 ;i<j ;i++){
            complexD prod = s1[i] * conj(s1[j]);
            P += (prod+conj(prod)) * sinc[i][j];  //plus rapide de * OU de conj()????????,
        }

        P += norm(s1[j]);
    }*/

    complexD P = 0;

    matrixC SVmulti = matrixC(fNumNus, vectorC(fNumNus, 0));
    for(int i = 0 ; i<fNumNus ; i++) {
        for(int j = 0 ; j<fNumNus ; j++){
            for(int k = 0 ; k<fNumNus ; k++){
                SVmulti[i][j] += fevolutionMatrixS[i][k] *V[k][j];
            }
        }
    }

    cout<<endl<<endl<<"E = "<<fEnergy<<endl;

    complexD sinc[fNumNus][fNumNus];
    for(int j = 0 ; j<fNumNus ; j++){
        for(int i = 0 ; i<j ; i++){
            double arg = (lambda[i] - lambda[j]) * dbin / 2; 
            sinc[i][j] = sin(arg) / arg;
            sinc[j][i] = sinc[i][j];

            cout<<"["<<i<<"]"<<"["<<j<<"] : "<<sinc[i][j]<<"    ";
        }

        sinc[j][j] = 1;
    }
    
    for(int j = 0 ; j<fNumNus ; j++){
        for(int i = 0 ;i<fNumNus ;i++){
            P += SVmulti[flvf][i] * conj(SVmulti[flvf][j]) * conj(V[flvi][i]) * V[flvi][j] * sinc[j][i];
        }
    }

    return real(P); 
}


//.............................................................................
///
///
///
double PMNS_TaylorExp::avgProbTaylor(int flvi, int flvf, double E , double dE)
{
    // reset K et S et Ve et lambdaE
    InitializeTaylorsVectors();

    SetEnergy(E);
    SetwidthBin(dE,0);

    //Propagate -> get S and K matrix (on the whole path)
    PropagateTaylor();

    //DiagolK -> get VE and lambdaE
    SolveK(fKE,flambdaE,fVE);

    //return fct avr proba
    return avgFormula(flvi,flvf,fdE*kGeV2eV, flambdaE, fVE);
}




//.............................................................................
///
/// 
///
double PMNS_TaylorExp::avgFormulaExtrapolation(int flvi, int flvf, double dbin, vectorD lambda, matrixC V)
{

    complexD P = 0;

    matrixC SVmulti = matrixC(fNumNus, vectorC(fNumNus, 0));
    for(int i = 0 ; i<fNumNus ; i++) {
        for(int j = 0 ; j<fNumNus ; j++){
            for(int k = 0 ; k<fNumNus ; k++){
                SVmulti[i][j] += fevolutionMatrixS[i][k] *V[k][j]; 
            }
        }
    }

    complexD exp[fNumNus][fNumNus];
    for(int i = 0 ; i<fNumNus ; i++){
        for(int j = 0 ; j<fNumNus; j++){
            double arg = (lambda[j] - lambda[i]) * dbin ; //DIVISER PAR2???    
            exp[j][i] = complexD(cos(arg),sin(arg));
        }
    }
    
    for(int j = 0 ; j<fNumNus ; j++){
        for(int i = 0 ;i<fNumNus ;i++){
            P += SVmulti[flvf][i] * conj(SVmulti[flvf][j]) * conj(V[flvi][i]) * V[flvi][j] * exp[j][i];
        }
    }

    //cout<<"P = "<<P<<endl;

    return real(P); 
}

//.............................................................................
///
///
///
double PMNS_TaylorExp::interpolationEnergy(int flvi, int flvf, double E , double dE)
{
    // reset K et S et Ve et lambdaE
    InitializeTaylorsVectors();

    SetEnergy(E);
    SetwidthBin(dE,0);

    //Propagate -> get S and K matrix (on the whole path)
    PropagateTaylor();

    //DiagolK -> get VE and lambdaE
    SolveK(fKE,flambdaE,fVE);

    return avgFormulaExtrapolation(flvi,flvf,fdE*kGeV2eV, flambdaE, fVE);
}
















//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
























//.............................................................................
///
///
///
vectorD PMNS_TaylorExp::ConvertLoEtoE(double LoE, double dLoE)
{
    // Make sure fPath is set
    // Use average if multiple paths
    SetCurPath(AvgPath(fNuPaths));

    vectorD Ebin(2);

    // Set a minimum energy
    double minLoE = 0.1 * LoE;

    // Transform range to E
    // Full range if low edge > minLoE
    if(LoE - dLoE / 2 >minLoE) {
        Ebin[0] = fPath.length * 0.5 * (1 / (LoE - dLoE / 2) + 1 / (LoE + dLoE / 2));
        Ebin[1] = fPath.length * (1 / (LoE - dLoE / 2) - 1 / (LoE + dLoE / 2));
    }
    else {
        Ebin[0] = fPath.length * 0.5 * (1 / minLoE + 1 / (LoE + dLoE / 2));
        Ebin[1] = fPath.length * (1 / minLoE - 1 / (LoE + dLoE / 2));
    }

    return Ebin;
}


//.............................................................................
///
///
///
double PMNS_TaylorExp::avgProbTaylorLoE(int flvi, int flvf, double LoE , double dLoE)
{
    if (LoE <= 0) return 0;

    if (fNuPaths.empty()) return 0;

    // Don't average zero width
    if (dLoE <= 0) return Prob(flvi, flvf, fPath.length / LoE);

    vectorD Ebin = ConvertLoEtoE(LoE,dLoE);

    return avgProbTaylor(flvi,flvf,Ebin[0],Ebin[1]);
}

//.............................................................................
///
///
///
double PMNS_TaylorExp::avgProbTaylorAngle(int flvi, int flvf, double E, double cosT , double dcosT)
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
    return avgFormula(flvi,flvf,fdcosT, flambdaCosT, fVcosT);
}

//.............................................................................
///
///
///
double PMNS_TaylorExp::avgAlgorithm(int flvi, int flvf)
{
    matrixC densityMatrix = matrixC(fNumNus, vectorC(fNumNus, 0));
    densityMatrix[flvi][flvi] = 1;

    RotateDensityM(true,fVcosT,densityMatrix);
    HadamardProduct(flambdaCosT,densityMatrix,fdcosT);
    RotateDensityM(false,fVcosT,densityMatrix);

    RotateDensityM(true,fVcosT,densityMatrix);
    HadamardProduct(flambdaE,densityMatrix,fdE);
    RotateDensityM(false,fVcosT,densityMatrix);

    RotateDensityM(false,fevolutionMatrixS,densityMatrix);

    return real(densityMatrix[flvf][flvf]);
}

//.............................................................................
///
///
///
double PMNS_TaylorExp::avgProbTaylor(int flvi, int flvf, double E , double dE, double cosT , double dcosT)
{
    // reset K et S et Ve et lambdaE
    InitializeTaylorsVectors();

    SetEnergy(E);
    SetCosT(cosT);
    SetwidthBin(dE,dcosT);

    //Propagate -> get S and K matrix (on the whole path)
    PropagateTaylor();

    //DiagolK -> get VE and lambdaE
    SolveK(fKE,flambdaE,fVE);
    SolveK(fKcosT,flambdaCosT,fVcosT);

    return avgAlgorithm(flvi,flvf);
}


