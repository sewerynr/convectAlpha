#ifndef FUN_INICJALIZ_H
#define FUN_INICJALIZ_H

#include "fvCFD.H"
//#include "fvOptions.H"

scalar funInitT1( double PsiZero, const dimensionedScalar& epsH, double par, const double x,const double y, const double z );

scalar funInitT4( double PsiZero, const dimensionedScalar& epsH, double par, const double x,const double y, const double z );

scalar funInitPsi1( double x, double y, double z );

scalar funInitPsi2( double x, double y, double z );

vector funInitV1( double x, double y, double z );

vector funInitV2( double x, double y, double z );

vector funInitV3( double x, double y, double z );

#endif // FUN_INICJALIZ_H
