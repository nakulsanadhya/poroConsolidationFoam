// wEqn.H — mechanics step: quasi-static equilibrium under uniaxial strain
//
//   d/dz( Mc dw/dz - alpha p ) = 0
{
    fvScalarMatrix wEqn
    (
        fvm::laplacian(Mc, w) == alpha*fvc::grad(p)().component(vector::Z)
    );

    wEqn.solve();
}

// ************************************************************************* //
