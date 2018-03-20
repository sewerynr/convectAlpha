#include "fun_initialize.h"

scalar funInitT1(double PsiZero, const dimensionedScalar& epsH, double par, const double x,const double y, const double z )
{
    scalar ret = 0.0;
    double zmieniacz = 4.0;
    //ret = 0.5*(scalar(1.) + Foam::tanh( PsiZero/(2.0*epsH.value()) ));
    ret = 1.0/(1.0+Foam::exp(-PsiZero/epsH.value()*zmieniacz));
    return ret;
}

/*
 * alpha(psi) circular interface
 */
scalar funInitT4( double PsiZero, const dimensionedScalar& epsH, double par, const double x,const double y, const double z )
{
    const static scalar R = 0.15;
    scalar eps = epsH.value();
    scalar ret = 0.0;
    scalar r = Foam::sqrt(x*x + (y-0.25)*(y-0.25));
    return 1 - 0.5*(1.0 + Foam::tanh((r-R)/2/eps));
}

scalar funInitPsi1( double x, double y, double z )
{
    scalar ret = 0.0;
    ret = x;
    return ret;
}

scalar funInitPsi2( double x, double y, double z )
{
    scalar ret = 0.0;
    scalar xc  = 0.0;
    scalar yc  = 0.0;
    ret = Foam::sqrt((x-xc)*x + (y-yc)*(y-yc));
    return ret;
}

vector funInitV1( double x, double y, double z )
{
    vector v(y*10.0,- x*10.0 , 0.0);
    return v;
}

vector funInitV2( double x, double y, double z )
{
//    using Foam::constant::mathematical::pi;
    using Foam::mathematicalConstant::pi;
    using Foam::sin;
    using Foam::cos;
    using Foam::pow;
    const static double v0 = 2.0; // powinno byc 1
    const static double L = 1.0;  // dl domeny
    const static double k = pi/L;

    x += 0.5;
    y += 0.5;

    return vector(-v0*sin(k*x)*sin(k*x)*sin(2.0*k*y), v0*sin(k*y)*sin(k*y)*sin(2.0*k*x), 0.0);
}
// odwrucony wektor predkosci
vector funInitV3( double x, double y, double z )
{
    using Foam::mathematicalConstant::pi;
    using Foam::sin;
    using Foam::cos;
    using Foam::pow;
    const static double v0 = 2.0;
    const static double L = 1.0; // dl domeny
    const static double k = pi/L;

    x += 0.5;
    y += 0.5;

    return vector(v0*sin(k*x)*sin(k*x)*sin(2.0*k*y), -v0*sin(k*y)*sin(k*y)*sin(2.0*k*x), 0.0);
}
