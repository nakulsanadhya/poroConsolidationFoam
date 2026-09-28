// Evaluate the coupled equations WITHOUT the fixed-stress iteration term.
// Mass is normalized as per-step fluid content / constrained storage*p0.
// Mechanics is force/volume normalized by sigma0/H (not a linear-solver residual).
{
    const volScalarField massDefect
    (
        "massDefect",
        storage*(p-pOld) + alpha*(eps-epsOld)
      - dt*fvc::laplacian(mobility, p)
      - fvc::laplacian(tau, p) + fvc::laplacian(tau, pOld)
    );
    const volScalarField mechanicsDefect
    (
        "mechanicsDefect",
        fvc::laplacian(Mc, w)-alpha*fvc::grad(p)().component(vector::Z)
    );
    massResidual = gMax(mag(massDefect.internalField()))/storageScale;
    mechanicsResidual = gMax(mag(mechanicsDefect.internalField()))/(sigma0.value()/H);
    forAll(p.internalField(), i)
    {
        if (!std::isfinite(p[i]) || !std::isfinite(w[i])
         || !std::isfinite(massDefect[i]) || !std::isfinite(mechanicsDefect[i]))
        {
            FatalErrorIn(args.executable()) << "NONFINITE_SOLUTION" << exit(FatalError);
        }
    }
}
