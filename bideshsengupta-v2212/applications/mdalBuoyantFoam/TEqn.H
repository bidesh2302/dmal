{
    volScalarField kappaEff
    (
        "kappaEff",
        turbulence->nu()/Pr + turbulence->nut()/Prt
    );

    fvScalarMatrix TEqn
    (
        fvm::ddt(T)
      + fvm::div(phi, T)
      - fvm::laplacian(kappaEff, T)
    );
    
    TEqn.solve();

    //~ rhok = 1.0 - beta*(T - TRef);
    rhok = beta*(T - TRef);
}
