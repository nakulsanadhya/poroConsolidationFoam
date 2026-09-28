// Integrated balance measured from this invocation's starting state.
// Restart segments intentionally have their own zero of cumulative drainage.
{
    scalar drainageRate = 0;
    const surfaceScalarField gradP("gradP", fvc::snGrad(p));
    forAll(mesh.boundary(), patchi)
    {
        if (mesh.boundary()[patchi].type() == "empty") continue;
        const scalarField& area = mesh.boundary()[patchi].magSf();
        forAll(area, f)
        {
            drainageRate -= mobility.value()*gradP.boundaryField()[patchi][f]*area[f];
        }
    }
    cumulativeDrain += dt.value()*drainageRate;
    const volScalarField physicalChange
    (
        "physicalChange", storage*(p-pStart)+alpha*(eps-epsStart)
    );
    const volScalarField stabilizationChange
    (
        "stabilizationChange", -fvc::laplacian(tau,p)+fvc::laplacian(tau,pStart)
    );
    const scalar physical = gSum(mesh.V()*physicalChange.internalField())/contentScale;
    const scalar stabilizationContent = gSum(mesh.V()*stabilizationChange.internalField())/contentScale;
    const scalar drained = cumulativeDrain/contentScale;
    const scalar balance = physical+stabilizationContent+drained;
    if (!std::isfinite(balance))
    {
        FatalErrorIn(args.executable()) << "NONFINITE_BALANCE" << exit(FatalError);
    }
    Info<< "DIAGNOSTICS " << runTime.value() << " " << massResidual << " "
        << mechanicsResidual << " " << physical << " " << stabilizationContent
        << " " << drained << " " << balance << endl;
}
