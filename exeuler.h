#ifndef EXEULER_H
#define EXEULER_H

#include <stdlib.h>
#include "localFunctions.h"

#include "fvCFD.H"
//#include "fvOptions.H"

void exEuler(const volScalarField& C, Time& runTime, const fvMesh& mesh, dimensionedScalar dtau,const volScalarField& PsiZero, volScalarField& Psi,
           volScalarField& T, volScalarField& Told, const dimensionedScalar& epsH, const double& eps, const bool& limitFieldT, const int& ilekrcz, double gamma, const bool& mapFunLog, const double ilePkt, const double gradPsiLimit);


#endif // EXEULER_H
