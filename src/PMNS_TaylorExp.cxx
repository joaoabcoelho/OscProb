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

    SetwidthBin(0.1);
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
void PMNS_TaylorExp::InitializeTaylorsVectors()
{
    fKE = matrixC(fNumNus, vectorC(fNumNus, 0));

    flambdaE = vectorD(fNumNus, 0);

    fVE = matrixC(fNumNus, vectorC(fNumNus, 0));

    fevolutionMatrixS = matrixC(fNumNus, vectorC(fNumNus, 0));
    for(int i= 0 ; i<fevolutionMatrixS.size(); i++){
        fevolutionMatrixS[i][i] = 1;
    }


    //REEL OU COMPLEX??????????
}

//.............................................................................
///
///
///
void PMNS_TaylorExp::SetwidthBin(double widthBin)
{
    fwidthBin = widthBin;
}

//.............................................................................
///
///
///
void PMNS_TaylorExp::rotateS(vectorC fPhases,matrixC& S)
{
    for(int j = 0 ; j<fNumNus ; j++) //PEUT ETRE AMELIORER
    {
        for(int i = 0 ; i<=j ; i++)
        {
            for(int k = 0 ; k<fNumNus ; k++)
            {
                S[i][j] += fEvec[i][k] * fPhases[k] * conj(fEvec[j][k]);
            }

            if(i != j){
                S[j][i] = -conj(S[i][j]);//??????????????????????
            }
        }
    }
}

//.............................................................................
///
///
///
void PMNS_TaylorExp::rotateK(matrixC Kmass,matrixC& Kflavor)
{
    for(int j = 0 ; j<fNumNus ; j++)
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
    }
}

//.............................................................................
///
///
///
void PMNS_TaylorExp::BuildKE(double L , matrixC& K)
{
    double lv = 2 * kGeV2eV * fEnergy; // 2E in eV

    for(int j = 0 ; j<fNumNus ; j++)
    {
        for(int i = 0 ; i<=j ; i++)
        {
            for(int k = 0 ; k<fNumNus ; k++)
            {
                for(int l = 0 ; l<fNumNus ; l++)
                {
                    K[i][j] += conj(fEvec[k][i]) * fHam[k][l] * fEvec[l][j];
                }
            }

            complexD C;
            if(i == j){
                C = complexD(1,0);
            }
            else{
                double argg = (fEval[i] - fEval[j]) * L;
                C = (complexD(cos(argg), sin(argg)) - complexD(1,0) ) / (complexD(0,argg)); 
                cout<<C<<endl;
            }// C=(1,0) because of H H' commutation (due to cst density case)

            K[i][j] *= (-L / lv) * K[i][j] * C; 

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
void PMNS_TaylorExp::MultiplicationRule(matrixC SLayer,matrixC KLayer)
{
    //K (ne pas modifier S avant de l'appliquer ici )
    matrixC KCopy = matrixC(fNumNus, vectorC(fNumNus, 0));
    for(int j = 0 ; j<fNumNus ; j++){
        for(int i = 0 ; i<fNumNus ; i++){ //PAS OBLIGIER DE COPIER
            KCopy[i][j] = fKE[i][j];
        }
    }
    for(int j = 0 ; j<fNumNus ; j++){
        for(int i = 0 ; i<=j ; i++){
            for(int k = 0 ; k<fNumNus ; k++){
                for(int l = 0 ; l<fNumNus ; l++){                                                          
                    fKE[i][j] += conj(fevolutionMatrixS[k][i]) * KLayer[k][l] * fevolutionMatrixS[l][j] + KCopy[i][j]; 
                }
            }

            if(i != j){
                fKE[j][i] = -conj(fKE[i][j]);
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

    matrixC SCopy = matrixC(fNumNus, vectorC(fNumNus, 0));
    for(int j = 0 ; j<fNumNus ; j++){
        for(int i = 0 ; i<fNumNus ; i++){
            SCopy[i][j] = fevolutionMatrixS[i][j];
        }
    }

    for(int j = 0 ; j<fNumNus ; j++){
        for(int i = 0 ; i<=j ; i++){
            for(int k = 0 ; k<fNumNus ; k++){

                fevolutionMatrixS[i][j] += SLayer[i][k] * SCopy[k][j];
            }

            if(i != j){
                fevolutionMatrixS[j][i] = -conj(fevolutionMatrixS[i][j]);
            }
        }

    }
    

}

//.............................................................................
///
///
///
double PMNS_TaylorExp::avrProbTaylor(int flvi, int flvf, double E , double widthBin)
{
    // reset K et S et Ve et lambdaE
    InitializeTaylorsVectors();

    SetEnergy(E);
    SetwidthBin(widthBin);

    //Propagate -> get S and K matrix (on the whole path)
    PropagateTaylor();

    //DiagolK -> get VE and lambdaE
    SolveK();

    //return fct avr proba
    return avgFormula(flvi,flvf);
}

//.............................................................................
///
/// Propagate neutrino state through full path
///
void PMNS_TaylorExp::PropagateTaylor()
{
  for (int i = 0; i < int(fNuPaths.size()); i++) { cout<<endl<<"path :"<<i<<endl;  PropagatePathTaylor(fNuPaths[i]); }
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
    double L = p.length;

    // Solve for eigensystem
    SolveHam();

    // Get the evolution matrix in mass basis
    double LengthIneV = kKm2eV * p.length;
    for (int i = 0; i < fNumNus; i++) {
        double arg = fEval[i] * LengthIneV;
        fPhases[i] = complexD(cos(arg), -sin(arg));
    }

    // Rotate S in flavor basis
    matrixC Sflavor = matrixC(fNumNus, vectorC(fNumNus, 0));
    rotateS(fPhases,Sflavor);                                   //probleme approx diag=1 

    // Build KE in mass basis
    matrixC Kmass = matrixC(fNumNus, vectorC(fNumNus, 0));
    BuildKE(L,Kmass);

    // Rotate KE in flavor basis
    matrixC Kflavor = matrixC(fNumNus, vectorC(fNumNus, 0));
    rotateK(Kmass,Kflavor);

    //multiplication rule for K and S 
    MultiplicationRule(Sflavor,Kflavor);



    for(int j=0 ; j<3 ; j++)
    {
        for(int k=0 ; k<3 ; k++)
        {
            cout<<fKE[j][k]<<" ";
        }
        cout<<endl;
    }
 
}

//.............................................................................
///
/// 
///
void PMNS_TaylorExp::SolveK()
{
    double   fEvalGLoBES[3];
    complexD fEvecGLoBES[3][3];

    complexD fKEconv[3][3];
    for(int j = 0 ; j<fNumNus ; j++){
        for(int i = 0 ; i<fNumNus ; i++){       //AAAAAAAAAAAAAAAAAA
            fKEconv[i][j] = fKE[i][j];
        }
    }

    // Solve Hamiltonian for eigensystem using the GLoBES method
    zheevh3(fKEconv, fEvecGLoBES, fEvalGLoBES);

    // Fill fEval and fEvec vectors from GLoBES arrays
    for (int i = 0; i < fNumNus; i++) {
        flambdaE[i] = fEvalGLoBES[i];
        for (int j = 0; j < fNumNus; j++) { fVE[i][j] = fEvecGLoBES[i][j]; }
    }
}

//.............................................................................
///
/// 
///
double PMNS_TaylorExp::avgFormula(int flvi, int flvf)
{
    vectorC SVmulti = vectorC(fNumNus, 0);
    for(int i = 0 ; i<fNumNus ; i++) {
        for(int j = 0 ; j<fNumNus ; j++){
            SVmulti[i] += fevolutionMatrixS[flvf][j] *fVE[j][i];
        }
    }

    complexD P;
    complexD s1[fNumNus]; //le plus rapide? vector+pushback OU tableau+size definie des le débuts?
    //complexD s2[fNumNus];
    complexD sinc[fNumNus][fNumNus];

    for(int i = 0 ; i<fNumNus ; i++){
        s1[i] = SVmulti[i] * conj(fVE[flvi][i]) ;
        //s2[i] = conj(SVmulti[i]) * fVE[flvi][i] ;  S1+CONJ(S2)!!!!!!!!!!!!!

        //sinc[i][i] = 1;
        for(int j = 0 ; j<i ; j++){
            double arg = (flambdaE[j] - flambdaE[i]) * fwidthBin;
            sinc[j][i] = sin(arg)/arg;
        }
    }

    for(int j = 0 ; j<fNumNus ; j++){
        for(int i = 0 ;i<j ;i++){
            complexD prod = s1[i] * conj(s1[j]);
            P += (prod+conj(prod)) * sinc[i][j];  //plus rapide de * OU de conj()????????,
        }

        P += norm(s1[j]);
    }

    cout<<"Proba = "<<P<<endl;

    return real(P); 
}