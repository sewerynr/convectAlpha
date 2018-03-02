# include "implicit.h"

void implicit(const volScalarField& C, Time& runTime, const fvMesh& mesh, dimensionedScalar dtau,const volScalarField& PsiZero, volScalarField& Psi,
           volScalarField& T, volScalarField& Told, const dimensionedScalar& epsH, const double& eps, const bool& limitFieldT, const int& ilkrcz, double gamma, const bool& mapFunLog)
{
    Info<< "Implicit !!! "  << endl;

    scalar one(1);
    dimensionedScalar SMALL_NUMBER("small", dimless, SMALL);

    Psi = epsH*Foam::log((T+eps)/(1-T+eps));
    surfaceScalarField phiR = createPhiFieldImp(C, runTime,mesh,T,PsiZero,Psi);

    fvScalarMatrix InitTEq
    (
        fvm::ddt(T)
        -fvm::div(phiR*dtau, T)
    );

    scalar res = Foam::max(Foam::mag(InitTEq.residual()));
    Info << "Start residual: " << res <<endl;
    for (int i = 1; i <= ilkrcz; ++i)
    {
        Told = T;
        fvScalarMatrix TEqn
        (
            fvm::ddt(T)
            -fvm::div(phiR*dtau, T)
        );

        TEqn.solve();

        if (limitFieldT)
            limitT(T);

        double norm1 = 0.;
        ++runTime;
        norm1 = Foam::sum(Foam::mag(T-Told)).value() / T.size();
        Info<< "Norma 1 = " << norm1 << endl;

        runTime.write();
        phiR = linearInterpolate( C*( scalar(1.) - T )*( mag(fvc::grad(Psi))- scalar(1.) ) * (fvc::grad(PsiZero) /( mag(fvc::grad(PsiZero)) + SMALL_NUMBER ))  )  & mesh.Sf();
        Psi = epsH*Foam::log((T+eps)/(1-T+eps));
    }
}
