// pEqn.H — flow step of the fixed-stress split (implicit Euler)
//
//   (1/M)(p - p^n)/dt + beta (p - p^k)/dt + alpha (eps^k - eps^n)/dt
//   - div(mobility grad p) - div(tau grad (p - p^n))/dt = 0
//
// eps^k is the strain from the latest mechanics solve; beta (p - p^k)/dt
// is the fixed-stress term and vanishes at convergence.
{
    fvScalarMatrix pEqn
    (
        fvm::Sp((storage + beta)/dt, p)
      - (storage/dt)*pOld
      - (beta/dt)*pIter
      + alpha*(eps - epsOld)/dt
      - fvm::laplacian(mobility, p)
      - fvm::laplacian(tau/dt, p)
      + fvc::laplacian(tau/dt, pOld)
    );

    pEqn.solve();
}

// ************************************************************************* //
