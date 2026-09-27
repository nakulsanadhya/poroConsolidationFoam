/*---------------------------------------------------------------------------*\
    poroConsolidationFoam

    A custom finite-volume solver for 1D Terzaghi consolidation:

        du/dt = cv * d2u/dz2

    where u is excess pore pressure and cv is the coefficient of
    consolidation (cv = k / (mv * gammaW)).

    Phase 1 (this file): the uncoupled consolidation (pressure diffusion)
    equation, verified against Terzaghi's exact series solution with
    mesh- and time-step-refinement studies (see verify.py). This checks
    the discretization and boundary conditions before adding true
    two-way coupling.

    Phase 2 (planned): couple this pressure equation to a displacement
    field via a fixed-stress split iteration (Kim, Tchelepi & Juanes,
    2011, "Stability and convergence of sequential methods for coupled
    flow and geomechanics"), giving genuine two-way poroelastic coupling.

    Independently written; not derived from any group/lab codebase.
    Built on foam-extend (open-source, GPL) as a standalone application.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << endl;

        // Implicit Euler in time, standard finite-volume Laplacian in space.
        // fvm::ddt / fvm::laplacian assemble the linear system; solve()
        // dispatches to the solver/tolerance set in fvSolution.
        fvScalarMatrix uEqn
        (
            fvm::ddt(u)
          - fvm::laplacian(cv, u)
        );

        uEqn.solve();

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
