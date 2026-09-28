/*---------------------------------------------------------------------------*\
    biotConsolidationFoam

    Two-way coupled 1D Biot consolidation (uniaxial strain), solved with the
    fixed-stress sequential split.

    Unknowns: pore pressure p [Pa] and vertical displacement w [m], both
    cell-centred. z points down, w is positive downward, stress is
    tension-positive.

      Equilibrium : d/dz( Mc dw/dz - alpha p ) = 0
      Fluid mass  : (1/M) dp/dt + alpha d(eps)/dt = d/dz( mobility dp/dz )
                    with eps = dw/dz

    Each time step iterates:
      1. flow step with the total stress frozen (fixed-stress term
         beta (p - p^k)/dt, beta = fixedStressFactor * alpha^2/Mc),
         using the volumetric strain from the latest mechanics solve;
      2. mechanics step for w with the new pressure;
    until the relative pressure change falls below fixedStressTolerance.

    Optional stabilization adds -d/dt div( tau grad p ) to the flow
    equation, tau = alpha^2 h^2 / (4 Mc). On this collocated grid the
    discrete strain equals (alpha/Mc)(p + h^2/4 Lap_h p) - sigma0/Mc, and the
    h^2/4 term makes the coupling blind to odd-even pressure modes; tau
    cancels it exactly (cf. Aguilar, Gaspar, Lisbona & Rodrigo, 2008).

    Independently written; not derived from any group/lab codebase.
    Built on foam-extend 4.1 (open source, GPL).
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include <cmath>
#include "fixedGradientFvPatchFields.H"

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    if (undrainedInitialState && mag(runTime.value()) < SMALL)
    {
        Info<< "Undrained initial pressure p0 = " << p0.value() << " Pa" << nl << endl;
        // Assignment leaves fixed-value patches (the drained top) unchanged
        p = p0;
        p.correctBoundaryConditions();
    }

    // Preserve saved displacement on restart; initialize mechanics only at t=0.
    if (mag(runTime.value()) < SMALL)
    {
        #include "wEqn.H"
    }
    eps = fvc::grad(w)().component(vector::Z);

    Info<< "Fixed-stress parameter beta = " << beta.value() << " 1/Pa" << nl
        << "Stabilization tau = " << tau.value() << " m^2/Pa" << nl << endl;

    const volScalarField pStart("pStart", p);
    const volScalarField epsStart("epsStart", eps);
    const scalar storageScale = (storage.value()+sqr(alpha.value())/Mc.value())*p0.value();
    const scalar contentScale = storageScale*totalVolume;
    scalar cumulativeDrain = 0;
    Info().precision(14);
    Info<< "Reliability checks v1: residuals and segment mass balance enabled" << endl;
    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        const dimensionedScalar dt = runTime.deltaT();
        const volScalarField pOld("pOld", p);
        const volScalarField epsOld("epsOld", eps);

        label iter = 0;
        scalar change = GREAT;
        scalar massResidual = GREAT;
        scalar mechanicsResidual = GREAT;
        if (!(dt.value() > 0) || !std::isfinite(dt.value()))
        {
            FatalErrorIn(args.executable()) << "Invalid time step" << exit(FatalError);
        }

        while (iter < maxFixedStressIterations)
        {
            iter++;
            const volScalarField pIter("pIter", p);

            #include "pEqn.H"
            #include "wEqn.H"
            eps = fvc::grad(w)().component(vector::Z);

            change = max(mag(p - pIter)).value()/p0.value();

            #include "equationResiduals.H"
            if (!std::isfinite(change) || !std::isfinite(massResidual)
             || !std::isfinite(mechanicsResidual))
            {
                FatalErrorIn(args.executable()) << "NONFINITE_SOLUTION" << exit(FatalError);
            }
            if (change < fixedStressTolerance && massResidual < equationTolerance
             && mechanicsResidual < equationTolerance)
            {
                break;
            }
        }

        Info<< "Fixed-stress iterations = " << iter
            << ", max relative pressure change = " << change << endl;

        if (change >= fixedStressTolerance || massResidual >= equationTolerance
         || mechanicsResidual >= equationTolerance)
        {
            FatalErrorIn(args.executable()) << "COUPLING_NOT_CONVERGED after " << iter
                << " iterations; change=" << change << " mass=" << massResidual
                << " mechanics=" << mechanicsResidual << exit(FatalError);
        }

        #include "massBalance.H"
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
