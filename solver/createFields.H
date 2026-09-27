// createFields.H
//
// Declares the excess pore pressure field u and reads the coefficient
// of consolidation cv from a transportProperties dictionary
// (constant/transportProperties), the same pattern OpenFOAM's own
// scalarTransportFoam/laplacianFoam tutorials use — kept deliberately
// familiar rather than reinventing conventions.

Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

// Coefficient of consolidation cv = k / (mv * gammaW)  [m^2/s]
dimensionedScalar cv
(
    transportProperties.lookup("cv")
);

Info<< "Reading field u (excess pore pressure)\n" << endl;

volScalarField u
(
    IOobject
    (
        "u",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// ************************************************************************* //
