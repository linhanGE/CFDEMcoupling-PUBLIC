/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright 2009-2012 JKU Linz
                                Copyright 2012-     DCS Computing GmbH, Linz
                                Copyright (C) 2013-     Graz University of  
                                                        Technology, IPPT
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling.

    CFDEMcoupling is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    This code is designed to realize coupled CFD-DEM simulations using LIGGGHTS
    and OpenFOAM(R). Note: this code is not part of OpenFOAM(R) (see DISCLAIMER).
\*---------------------------------------------------------------------------*/

#include "error.H"

#include "constDiffCentreSmoothing.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(constDiffCentreSmoothing, 0);

addToRunTimeSelectionTable
(
    smoothingModel,
    constDiffCentreSmoothing,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
constDiffCentreSmoothing::constDiffCentreSmoothing
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    smoothingModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    minAlphas_(readScalar(propsDict_.lookup("minAlphas"))),
    maxAlphas_(readScalar(propsDict_.lookup("maxAlphas"))),
    smoothingLength_(dimensionedScalar("smoothingLength",dimensionSet(0,1,0,0,0,0,0), readScalar(propsDict_.lookup("smoothingLength")))),
    DT_("DT", dimensionSet(0,2,-1,0,0), 0.),
    verbose_(false)
{

    if(propsDict_.found("verbose"))  
        verbose_ = true;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

constDiffCentreSmoothing::~constDiffCentreSmoothing()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
bool constDiffCentreSmoothing::doSmoothing() const
{
    return true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void Foam::constDiffCentreSmoothing::smoothen(volScalarField& fieldSrc) const
{
	volScalarField diffWorkField
    (
        IOobject
        (
            "tempDiffScalar",
            particleCloud_.mesh().time().timeName(),
            particleCloud_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        particleCloud_.mesh(),
        dimensionedScalar("zero", dimensionSet(0,0,0,0,0), 0),
        zeroGradientFvPatchScalarField::typeName
    );
    // Create scalar smooth field from virgin scalar smooth field template
    diffWorkField.primitiveFieldRef() = fieldSrc.primitiveFieldRef();

    double deltaT = diffWorkField.mesh().time().deltaTValue();
    DT_.value() = smoothingLength_.value() * smoothingLength_.value() / deltaT;
    
	Info<< "smoothing " << fieldSrc.name() << endl;
    // do smoothing
    solve
    (
        fvm::ddt(diffWorkField)
       -fvm::laplacian(DT_, diffWorkField)
    );

    // bound sSmoothField_
    forAll(diffWorkField.primitiveFieldRef(),cellI)
    {
        if (diffWorkField.primitiveFieldRef()[cellI] > maxAlphas_)
			Info << "Unphysical alphas found" << endl;		
		diffWorkField.primitiveFieldRef()[cellI]=max(minAlphas_,min(maxAlphas_,diffWorkField.primitiveFieldRef()[cellI]));
    }  

    // get data from working sSmoothField - will copy only values at new time
    fieldSrc.primitiveFieldRef() = diffWorkField.primitiveFieldRef();
    fieldSrc.correctBoundaryConditions(); 

    if(verbose_)
    {
        Info << "min/max(fieldSrc): " << min(fieldSrc) << tab << max(fieldSrc) << endl;
    }

}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void Foam::constDiffCentreSmoothing::smoothen(volVectorField& fieldSrc) const
{
    volVectorField diffWorkField
    (
        IOobject
        (
            "tempDiffVector",
            particleCloud_.mesh().time().timeName(),
            particleCloud_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        particleCloud_.mesh(),
        dimensionedVector("zero", dimensionSet(0,0,0,0,0), vector::zero),
        zeroGradientFvPatchVectorField::typeName
    );
	// Create scalar smooth field from virgin scalar smooth field template
    diffWorkField.primitiveFieldRef() = fieldSrc.primitiveFieldRef(); 

    double deltaT = diffWorkField.mesh().time().deltaTValue();
    DT_.value() = smoothingLength_.value() * smoothingLength_.value() / deltaT;
    
	Info<< "smoothing " << fieldSrc.name() << endl;
    // do smoothing
    solve
    (
        fvm::ddt(diffWorkField)
       -fvm::laplacian(DT_, diffWorkField)
    );

    // get data from working vSmoothField
    fieldSrc.primitiveFieldRef() = diffWorkField.primitiveFieldRef();
    fieldSrc.correctBoundaryConditions();  

    if(verbose_)
    {
        Info << "min/max(fieldSrc): " << min(fieldSrc) << tab << max(fieldSrc) << endl;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void Foam::constDiffCentreSmoothing::UsSmoothen(volVectorField& fieldSrc, volScalarField& alphas) const
{
    // Create scalar smooth field from virgin scalar smooth field template
    volVectorField diffWorkField
    (
        IOobject
        (
            "tempDiffVector",
            particleCloud_.mesh().time().timeName(),
            particleCloud_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        particleCloud_.mesh(),
        dimensionedVector("zero", dimensionSet(0,0,0,0,0), vector::zero),
        zeroGradientFvPatchVectorField::typeName
    );
	
    diffWorkField.primitiveFieldRef() = fieldSrc.primitiveFieldRef(); 

    double deltaT = diffWorkField.mesh().time().deltaTValue();
    DT_.value() = smoothingLength_.value() * smoothingLength_.value() / deltaT;

    // do the smoothing
    solve
    (
        fvm::ddt(diffWorkField)
       -fvm::laplacian(DT_, diffWorkField)
    );
	
	// get data from working diffWorkField - will copy only values at new time
	forAll(diffWorkField.primitiveFieldRef(), cellI)
    {
        if (alphas.primitiveFieldRef()[cellI] > ROOTVSMALL)
        {
            diffWorkField.primitiveFieldRef()[cellI] /=
                (alphas.primitiveFieldRef()[cellI]);
        }
    }

    // get data from working vSmoothField
    fieldSrc.primitiveFieldRef() = diffWorkField.primitiveFieldRef();
    fieldSrc.correctBoundaryConditions(); 

    if(verbose_)
    {
        Info << "min/max(fieldSrc): " << min(fieldSrc) << tab << max(fieldSrc) << endl;
    }

}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::constDiffCentreSmoothing::smoothenReferenceField(volVectorField& fieldSrc) const   // for pure virtual function
{
	// do nothing
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
