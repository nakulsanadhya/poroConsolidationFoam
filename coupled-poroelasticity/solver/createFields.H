// createFields.H — biotConsolidationFoam
//
// Reads the poroelastic parameters and solution controls from
// constant/biotProperties, and the fields p (pore pressure) and
// w (vertical displacement).

Info<< "Reading biotProperties\n" << endl;

IOdictionary biotProperties
(
    IOobject
    (
        "biotProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

// Poroelastic parameters
const dimensionedScalar Mc(biotProperties.lookup("Mc"));            // drained constrained modulus K + 4G/3 [Pa]
const dimensionedScalar M(biotProperties.lookup("M"));              // Biot modulus [Pa]
const dimensionedScalar alpha(biotProperties.lookup("alpha"));      // Biot coefficient [-]
const dimensionedScalar mobility(biotProperties.lookup("mobility"));// k/mu [m^2/(Pa s)]
const dimensionedScalar sigma0(biotProperties.lookup("sigma0"));    // compressive load on loadedPatch [Pa]
const word loadedPatch(biotProperties.lookup("loadedPatch"));       // drained, loaded boundary

// Solution controls
const Switch undrainedInitialState(biotProperties.lookup("undrainedInitialState"));
const Switch stabilization(biotProperties.lookup("stabilization"));
const scalar fixedStressFactor(readScalar(biotProperties.lookup("fixedStressFactor")));
const scalar fixedStressTolerance(readScalar(biotProperties.lookup("fixedStressTolerance")));
const label maxFixedStressIterations(readLabel(biotProperties.lookup("maxFixedStressIterations")));

const scalar equationTolerance = biotProperties.found("equationTolerance")
    ? readScalar(biotProperties.lookup("equationTolerance")) : 1e-8;
if (!std::isfinite(Mc.value()) || !std::isfinite(M.value())
 || !std::isfinite(mobility.value()) || !std::isfinite(sigma0.value())
 || !std::isfinite(alpha.value()) || !std::isfinite(fixedStressFactor)
 || !std::isfinite(fixedStressTolerance) || !std::isfinite(equationTolerance)
 || !(Mc.value() > 0 && M.value() > 0 && mobility.value() > 0
   && sigma0.value() > 0 && alpha.value() > 0 && alpha.value() <= 1
   && fixedStressFactor > 0 && fixedStressTolerance > 0
   && equationTolerance > 0 && maxFixedStressIterations > 0))
{
    FatalErrorIn(args.executable()) << "Invalid material or solver controls"
        << exit(FatalError);
}

// Derived quantities
const dimensionedScalar storage = dimensionedScalar("one", dimless, 1.0)/M;  // 1/M [1/Pa]
const dimensionedScalar beta = fixedStressFactor*sqr(alpha)/Mc;     // fixed-stress term [1/Pa]
const dimensionedScalar p0 = alpha*M*sigma0/(Mc + sqr(alpha)*M);    // undrained pressure [Pa]

// The exact stabilization identity is restricted to this serial uniform column.
if (Pstream::parRun() || mesh.nCells() < 2 || mesh.nInternalFaces() != mesh.nCells()-1)
{
    FatalErrorIn(args.executable()) << "Requires a serial 1D column with >=2 cells"
        << exit(FatalError);
}
// Cell size for the stabilization term (checked below)
const dimensionedScalar h
(
    "h", dimLength, 1.0/gMax(mesh.deltaCoeffs().internalField())
);
const dimensionedScalar tau =
    (stabilization ? 1.0 : 0.0)*sqr(alpha)*sqr(h)/(4.0*Mc);         // [m^2/Pa]

Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field w\n" << endl;
volScalarField w
(
    IOobject
    (
        "w",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "validateCase.H"

// Traction on the loaded (drained) boundary: Mc dw/dz - alpha p = -sigma0
// with p = 0 there, i.e. a fixed normal gradient dw/dn = sigma0/Mc (the
// outward normal of the top boundary points in -z).
{
    const label patchi = mesh.boundaryMesh().findPatchID(loadedPatch);
    if (patchi < 0)
    {
        FatalErrorIn(args.executable())
            << "loadedPatch " << loadedPatch << " not found" << exit(FatalError);
    }
    fixedGradientFvPatchScalarField& wLoaded =
        refCast<fixedGradientFvPatchScalarField>(w.boundaryField()[patchi]);
    wLoaded.gradient() = sigma0.value()/Mc.value();
    w.correctBoundaryConditions();
}

// Volumetric strain eps = dw/dz (Gauss linear, same as the strain in the
// mass balance)
volScalarField eps
(
    IOobject
    (
        "eps",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::grad(w)().component(vector::Z)
);

// ************************************************************************* //
