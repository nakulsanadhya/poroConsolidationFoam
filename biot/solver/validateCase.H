// Deliberately restrictive domain: uniform, orthogonal, constant-area z column.
const scalar H = h.value()*mesh.nCells();
const scalar totalVolume = gSum(mesh.V());
{
    const scalar epsGeom = 1e-8;
    const vectorField& centres = mesh.C().internalField();
    const scalar cellV = totalVolume/mesh.nCells();
    forAll(centres, i)
    {
        if (mag(centres[i].x()-centres[0].x()) > epsGeom*H
         || mag(centres[i].y()-centres[0].y()) > epsGeom*H
         || mag(mesh.V()[i]/cellV-1) > epsGeom)
        {
            FatalErrorIn(args.executable()) << "Nonuniform/non-column mesh unsupported"
                << exit(FatalError);
        }
    }
    forAll(mesh.deltaCoeffs().internalField(), f)
    {
        const vector d = centres[mesh.neighbour()[f]]-centres[mesh.owner()[f]];
        const vector sf = mesh.Sf()[f];
        if (mag(mesh.deltaCoeffs()[f]*h.value()-1) > epsGeom
         || mag(mag(d.z())/h.value()-1) > epsGeom
         || mag(sf.x())+mag(sf.y()) > epsGeom*mag(sf)
         || mag(mag(sf)*h.value()/cellV-1) > epsGeom)
        {
            FatalErrorIn(args.executable()) << "Nonorthogonal/nonuniform mesh unsupported"
                << exit(FatalError);
        }
    }
    if (loadedPatch != "drained" || mesh.boundary().size() != 3)
    {
        FatalErrorIn(args.executable()) << "Expected drained, impermeable, sides patches"
            << exit(FatalError);
    }
    forAll(mesh.boundary(), patchi)
    {
        const word name = mesh.boundary()[patchi].name();
        const word pt = p.boundaryField()[patchi].type();
        const word wt = w.boundaryField()[patchi].type();
        if (name == "sides")
        {
            if (mesh.boundary()[patchi].type() != "empty" || pt != "empty" || wt != "empty")
            {
                FatalErrorIn(args.executable()) << "sides must be empty" << exit(FatalError);
            }
            const vectorField& sf = mesh.boundary()[patchi].Sf();
            forAll(sf, f)
            {
                if (mag(sf[f].z()) > epsGeom*mag(sf[f]))
                {
                    FatalErrorIn(args.executable()) << "Side faces must be parallel to z"
                        << exit(FatalError);
                }
            }
            continue;
        }
        const bool top = name == "drained";
        if ((!top && name != "impermeable") || mesh.boundary()[patchi].size() != 1
         || pt != (top ? "fixedValue" : "zeroGradient")
         || wt != (top ? "fixedGradient" : "fixedValue"))
        {
            FatalErrorIn(args.executable()) << "Unsupported boundary conditions"
                << exit(FatalError);
        }
        const vector sf = mesh.boundary()[patchi].Sf()[0];
        const scalar dc = mesh.boundary()[patchi].deltaCoeffs()[0];
        if (mag(sf.x())+mag(sf.y()) > epsGeom*mag(sf)
         || (top ? sf.z() >= 0 : sf.z() <= 0)
         || mag(dc*h.value()/2-1) > epsGeom
         || mag(mag(sf)*h.value()/cellV-1) > epsGeom
         || (top ? mag(p.boundaryField()[patchi][0]) > SMALL
                 : mag(w.boundaryField()[patchi][0]) > SMALL))
        {
            FatalErrorIn(args.executable()) << "Expected zero drainage pressure, fixed base, half-cell faces"
                << exit(FatalError);
        }
    }
}

// The cancellation proof uses these defaults, without field-specific overrides.
{
    const dictionary& schemes = mesh.schemesDict();
    const char* sections[] = {"gradSchemes", "laplacianSchemes", "interpolationSchemes", "snGradSchemes"};
    const char* expected[][3] = {{"Gauss", "linear", ""}, {"Gauss", "linear", "corrected"},
                                {"linear", "", ""}, {"corrected", "", ""}};
    const label lengths[] = {2,3,1,1};
    for (label i=0; i<4; ++i)
    {
        const dictionary& section = schemes.subDict(sections[i]);
        if (section.size() != 1 || !section.found("default"))
        {
            FatalErrorIn(args.executable()) << "Unsupported scheme overrides in " << sections[i]
                << exit(FatalError);
        }
        ITstream& entry = section.lookup("default");
        if (entry.size() != lengths[i])
        {
            FatalErrorIn(args.executable()) << "Unsupported scheme in " << sections[i]
                << exit(FatalError);
        }
        for (label j=0; j<lengths[i]; ++j)
        {
            const word actual(entry);
            if (actual != expected[i][j])
            {
                FatalErrorIn(args.executable()) << "Unsupported scheme in " << sections[i]
                    << exit(FatalError);
            }
        }
        entry.rewind();
    }
}
