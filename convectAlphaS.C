/*---------------------------------------------------------------------------*\
| OpenFOAM: The Open Source CFD Toolbox
\*---------------------------------------------------------------------------*/
#include "fvCFD.H"
#include "simpleControl.H"
#include <fstream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#include "localFunctions.h"
#include "exeuler.h"
#include "exrk3.h"
#include "exrk3face.h"
#include "implicit.h"
#include "fun_initialize.h"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#include "setRootCase.H"   // check case dir
#include "createTime.H"
#include "createMesh.H"
//    simpleControl simple(mesh);
#include "createFields.H"
//#include "createPhi.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#include "CourantNo.H"
#include "settings.h"

    dimensionedScalar dx( "dx", dimLength, Foam::pow(Foam::average(mesh.cellVolumes()), one/3.0));
    dimensionedScalar epsH( "epsH", dimLength, dx.value() / parEpsH);
    Info << "Ntau= " << Ntau << endl;
    epsH = epsH*Foam::sqrt(pDim);

    dimensionedScalar dtau( "dTau", dimTime, epsH.value() / pardTau );
    dimensionedScalar dtauI( "dTauI", dimless, epsH.value() / pardTau );

//  InitPsiXYZ(PsiZero, mesh, funInitPsi2); //2D
    InitPsiXYZ(PsiZero, mesh, funInitPsi1); //1D

//  compute T from PsiZero
    InitT(T, PsiZero, mesh, epsH, funInitT1, par);

    Info << endl << "current runTime=" << runTime.value() << endl;
    T.write();
    Info << "dx   =  "      << dx.value() << endl;
    Info << "epsH =  "      << epsH.value() << endl;
    Info << "dtau =  "      << dtau.value() << endl;
    Info << "ilosc Pkt =  " << ilePkt << endl << endl;

//    exRK3(C, runTime, mesh, dtau, PsiZero, Psi, T, Told, epsH, eps, limitFieldT, Ntau, gamma, ilePkt, gradPsiLimit);
//    exEuler(C, runTime, mesh, dtau, PsiZero, Psi, T, Told, epsH, eps, limitFieldT, Ntau, gamma, mapFunLog, ilePkt, gradPsiLimit);

//    exRK3Face(C, runTime, mesh, dtau, PsiZero, Psi, T, Told, epsH, eps, limitFieldT, Ntau, gamma, mapFunLog, ilePkt, gradPsiLimit);
    implicit(C, runTime, mesh, dtauI, PsiZero, Psi, T, Told, epsH, eps, limitFieldT, Ntau, gamma, ilePkt, gradPsiLimit);


    Info<< "End\n" << endl;
    return 0;
}
