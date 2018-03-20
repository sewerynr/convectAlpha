#ifndef CREATEFIELDSEXRK3_H
#define CREATEFIELDSEXRK3_H

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
    mesh,
    dimensionedScalar("gradPsi", dimless, scalar(-100000000.0))
);


surfaceScalarField msnGradPsi
(
    IOobject
    (
        "msngradPsi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("msngradPsi", dimless, scalar(-100000000.0))
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
    mesh,
    dimensionedScalar("gradTanalit", dimless/dimLength, scalar(-100000000.0))
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
    mesh,
    dimensionedScalar("magGradT", dimless/dimLength, scalar(-100000000.0))
);

volScalarField LaplaceTStraignt
(
    IOobject
    (
        "LaplaceTStraignt",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("LaplaceTStraignt", dimless/(dimLength*dimLength), scalar(-100000000.0))
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
    mesh,
    dimensionedScalar("TAnalit", dimless, scalar(-100000000.0))
    //0.5*(1.+Foam::tanh(PsiZero/(2.*epsH)))
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
    mesh,
    dimensionedScalar("LaplaceTW1D", dimless/(dimLength*dimLength), scalar(-100000000.0))
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
    mesh,
    dimensionedScalar("LaplaceAnalit", dimless/(dimLength*dimLength), scalar(-100000000.0))
 );

surfaceScalarField surfgradT
(   // nie grad interp lin na scianki tego co w nawiasie
    IOobject
    (
        "surfgradT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("surfgradT", (dimLength), scalar(-100000000.0))
);


volScalarField GradTW
(
    IOobject
    (
        "GradTW",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("surfgradT", dimless/(dimLength), scalar(-100000000.0))
);

surfaceScalarField phiR
(
    IOobject
    (
        "phiR",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("phiR", (dimLength*dimLength*dimLength)/dimTime, scalar(-100000000.0))
);

#endif // CREATEFIELDSEXRK3_H
