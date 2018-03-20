#ifndef CREATEFIELDSEXRK3FACE_H
#define CREATEFIELDSEXRK3FACE_H

surfaceScalarField phiR2
(
    IOobject
    (
        "phiR2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("phiR2", (dimLength*dimLength*dimLength)/dimTime, scalar(1.))
);

volScalarField mGradPsi
(
    IOobject
    (
        "gradPsi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mag(fvc::grad(Psi))
);

volScalarField gradTanalit
(
    IOobject
    (
        "gradTanalit",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mag( ( one - Foam::tanh(PsiZero/(2*epsH))*Foam::tanh(PsiZero/(2*epsH)) )*(1/(4*epsH)) )
);

volScalarField magGradT
(
    IOobject
    (
        "magGradT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mag(fvc::grad(T))
);

volVectorField GradT
(
    IOobject
    (
        "GradT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    T*(1-T)*(fvc::grad(Psi))/epsH
);

volScalarField sGradT
(
    IOobject
    (
        "sGradT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mag(fvc::grad(Psi))/epsH
);

volScalarField LaplaceT
(
    IOobject
    (
        "LaplaceT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    T*(1-T)*(fvc::laplacian(Psi))/epsH
);

volScalarField LaplaceTStraight
(
    IOobject
    (
        "LaplaceTStraight",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::laplacian(T)
);

volScalarField TAnalit
(
    IOobject
    (
        "TAnalit",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    0.5*(1.+Foam::tanh(PsiZero/(2.*epsH)))
);
volScalarField TStartowe
(
    IOobject
    (
        "TStartowe",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    T
);

volScalarField LaplaceTW1D
(
    IOobject
    (
        "LaplaceTW1D",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    T*(1.-T)*( fvc::laplacian(Psi) + (1./epsH)*( fvc::grad(Psi)&fvc::grad(Psi) )*(1.-2.*T) )/epsH
 );

volScalarField LaplaceAnalit
(
    IOobject
    (
        "LaplaceAnalit",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
        - ( Foam::tanh(PsiZero/(2.0*epsH)) - Foam::tanh(PsiZero/(2.0*epsH))*Foam::tanh(PsiZero/(2.0*epsH))*Foam::tanh(PsiZero/(2.0*epsH)) )/(4*epsH*epsH)
 );

surfaceScalarField surfgradT
(
    IOobject
    (
        "surfgradT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
   linearInterpolate( T*(1-T)*(fvc::grad(Psi))/epsH ) & mesh.Sf()
);

volVectorField GradTW
(
    IOobject
    (
        "GradTW",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    T*(1.-T)*(fvc::grad(Psi))/epsH
);

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   NEW  FOR FACE BASED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
surfaceScalarField PsiF
(
    IOobject
    (
        "PsiF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("PsiF", dimLength, scalar(0))
);

surfaceScalarField TF
(
    IOobject
    (
        "TF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("TF", dimless, scalar(0))
);

surfaceScalarField Cf
(
    IOobject
    (
        "Cf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("Cf", dimLength/dimTime, scalar(1.))
);

surfaceScalarField mGradPsiF
(
    IOobject
    (
        "mGradPsiF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::snGrad(Psi)
);

#endif // CREATEFIELDSEXRK3FACE_H
