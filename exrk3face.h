#ifndef EXRK3FACE_H
#define EXRK3FACE_H


#include <stdlib.h>
#include "localFunctions.h"

#include "fvCFD.H"

void exRK3Face(const volScalarField& C, Time& runTime, const fvMesh& mesh, dimensionedScalar dtau,const volScalarField& PsiZero, volScalarField& Psi,
           volScalarField& T, volScalarField& Told, const dimensionedScalar& epsH, const double& eps, const bool& limitFieldT, const int& ilekrcz, double gamma, const bool& mapFunLog, const double ilePkt, const double gradPsiLimit);


#endif // EXRK3FACE_H

