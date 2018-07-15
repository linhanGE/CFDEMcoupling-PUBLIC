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

#include "diffusionSmoothing.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(diffusionSmoothing, 0);

addToRunTimeSelectionTable
(
    smoothingModel,
    diffusionSmoothing,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
diffusionSmoothing::diffusionSmoothing
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    smoothingModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    minAlphas_(readScalar(propsDict_.lookup("minAlphas"))),
    maxAlphas_(readScalar(propsDict_.lookup("maxAlphas"))),
	diffusionBandWidth_(readScalar(propsDict_.lookup("diffusionBandWidth"))),
	diffusionSteps_(readLabel(propsDict_.lookup("diffusionSteps"))),
	smoothDirection_(propsDict_.lookup("smoothDirection")),
    verbose_(false),
    /******************************************/
    /* diffusion method */
	diffusionRunTime_
    (
        "controlDiffDict",
        sm.mesh().time().rootPath(),
        sm.mesh().time().caseName()
    ),
    diffusionMesh_
    (
        IOobject
        (
            fvMesh::defaultRegion,
            diffusionRunTime_.timeName(),
            diffusionRunTime_,
            IOobject::MUST_READ
        )
    ),
//******************************************//
	simple_(diffusionMesh_),
    diffusionTimeCount_(2,0.0)  // length is 2, initial value is 0
{
    if(propsDict_.found("verbose"))  
        verbose_ = true;
	// determine the time and time step in diffusion procedure
    scalar diffusionTime = pow(diffusionBandWidth_, 2)/4;
    scalar diffusionDeltaT = diffusionTime/(diffusionSteps_ + ROOTVSMALL);

    diffusionRunTime_.setEndTime(diffusionTime);
    diffusionRunTime_.setDeltaT(diffusionDeltaT);
    Info << "diffusion time is: " << diffusionTime << endl;
    Info << "diffusion time step is: " << diffusionDeltaT << endl;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

diffusionSmoothing::~diffusionSmoothing()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
bool diffusionSmoothing::doSmoothing() const
{
    return true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void Foam::diffusionSmoothing::smoothen(volScalarField& fieldSrc) const
{
	scalar t0 = particleCloud_.mesh().time().elapsedCpuTime();
	volScalarField diffWorkField
    (
        IOobject
        (
            "tempDiffScalar",
            diffusionRunTime_.timeName(),
            diffusionMesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        diffusionMesh_,
        dimensionedScalar
        (
            "zero",
            dimless,
            scalar(0.0)
        ),
        zeroGradientFvPatchScalarField::typeName
    );

	diffWorkField.primitiveFieldRef() = fieldSrc.primitiveFieldRef();
	
	dimensionedTensor DT("DT", dimensionSet(0, 2, -1, 0, 0), smoothDirection_);
	
	scalar startTime = diffusionRunTime_.startTimeIndex();
    label startIndex = diffusionRunTime_.timeIndex();

	Info<< "smoothing " << fieldSrc.name() << endl;
	
	diffusionTimeCount_[0] += particleCloud_.mesh().time().elapsedCpuTime() - t0;
    t0 = particleCloud_.mesh().time().elapsedCpuTime();
	
	while (diffusionRunTime_.loop())
    {
        if (diffusionRunTime_.timeIndex() == 1)
        {
            while (simple_.correctNonOrthogonal())
            {
                solve(fvm::ddt(diffWorkField) - fvm::laplacian(DT, diffWorkField));
            }
        }
        else
        {
            solve(fvm::ddt(diffWorkField) - fvm::laplacian(DT, diffWorkField));
        }
    }

    diffusionTimeCount_[1] += particleCloud_.mesh().time().elapsedCpuTime() - t0;
    t0 = particleCloud_.mesh().time().elapsedCpuTime();

    diffusionRunTime_.setTime(startTime,startIndex);
    
	// bound diffWorkField
    forAll(diffWorkField.primitiveFieldRef(),cellI)
    {
        diffWorkField.primitiveFieldRef()[cellI]=max(minAlphas_,min(maxAlphas_,diffWorkField.primitiveFieldRef()[cellI]));
    }
	
	// get data from working diffWorkField - will copy only values at new time
    fieldSrc.primitiveFieldRef() = diffWorkField.primitiveFieldRef();
    diffusionTimeCount_[0] += particleCloud_.mesh().time().elapsedCpuTime() - t0;
	fieldSrc.correctBoundaryConditions(); 

    if(verbose_)
    {

        Info << "min/max(alphas): " << min(fieldSrc) << tab << max(fieldSrc) << endl;
    }

}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::diffusionSmoothing::smoothen(volVectorField& fieldSrc) const
{
	scalar t0 = particleCloud_.mesh().time().elapsedCpuTime();
    volVectorField diffWorkField
    (
        IOobject
        (
            "tempDiffVector",
            diffusionRunTime_.timeName(),
            diffusionMesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        diffusionMesh_,
        dimensionedVector
        (
            "zero",
            dimless,
            vector::zero
        ),
        zeroGradientFvPatchVectorField::typeName
    );
	
    diffWorkField.primitiveFieldRef() = fieldSrc.primitiveFieldRef(); 
	
    dimensionedTensor DT("DT", dimensionSet(0, 2, -1, 0, 0), smoothDirection_);
	
	scalar startTime = diffusionRunTime_.startTimeIndex();
    label startIndex = diffusionRunTime_.timeIndex();

	Info<< "smoothing " << fieldSrc.name() << endl;
    vector Ftotal1(vector::zero);
    if(verbose_) 
    {
        forAll(fieldSrc.primitiveFieldRef(), cellI)
        {
            Ftotal1 += fieldSrc.primitiveFieldRef()[cellI];
        }
    }

	diffusionTimeCount_[0] += particleCloud_.mesh().time().elapsedCpuTime() - t0;
    t0 = particleCloud_.mesh().time().elapsedCpuTime();
	
	while (diffusionRunTime_.loop())
    {
        Info << "timeIndex = " << diffusionRunTime_.timeIndex() << endl;
        if (diffusionRunTime_.timeIndex() == 1)
        {
            while (simple_.correctNonOrthogonal())
            {
                solve(fvm::ddt(diffWorkField) - fvm::laplacian(DT, diffWorkField));
            }
        }
        else
        {
            solve(fvm::ddt(diffWorkField) - fvm::laplacian(DT, diffWorkField));
        }
    }

    diffusionTimeCount_[1] += particleCloud_.mesh().time().elapsedCpuTime() - t0;
    t0 = particleCloud_.mesh().time().elapsedCpuTime();

    diffusionRunTime_.setTime(startTime,startIndex);

    // get data from working diffWorkField - will copy only values at new time
    fieldSrc.primitiveFieldRef() = diffWorkField.primitiveFieldRef();
    vector Ftotal2(vector::zero);
    if(verbose_) 
    {
        forAll(fieldSrc.primitiveFieldRef(), cellI)
        {
            Ftotal2 += fieldSrc.primitiveFieldRef()[cellI];
        }
    }
    
    diffusionTimeCount_[0] += particleCloud_.mesh().time().elapsedCpuTime() - t0;

    fieldSrc.correctBoundaryConditions(); 

    if(verbose_)
    {
        reduce(Ftotal1,sumOp<vector>());
        reduce(Ftotal2,sumOp<vector>());
        Info << "total f before smoothing: " << Ftotal1 << endl;
        Info << "total f after smoothing: " << Ftotal2 << endl;
        Info << "min/max(f): " << min(fieldSrc) << tab << max(fieldSrc) << endl;
    }
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::diffusionSmoothing::UsSmoothen(volVectorField& fieldSrc, volScalarField& alphas) const
{
    scalar t0 = particleCloud_.mesh().time().elapsedCpuTime();
    volVectorField diffWorkField
    (
        IOobject
        (
            "tempDiffVector",
            diffusionRunTime_.timeName(),
            diffusionMesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        diffusionMesh_,
        dimensionedVector
        (
            "zero",
            dimless,
            vector::zero
        ),
        zeroGradientFvPatchVectorField::typeName
    );
	
    diffWorkField.primitiveFieldRef() = fieldSrc.primitiveFieldRef(); 
	
    dimensionedTensor DT("DT", dimensionSet(0, 2, -1, 0, 0), smoothDirection_);
	
	scalar startTime = diffusionRunTime_.startTimeIndex();  // Return start time index
    label startIndex = diffusionRunTime_.timeIndex();       // Return current time index. 

	Info<< "smoothing " << fieldSrc.name() << "*" << alphas.name() << endl;

    vector Utotal1(vector::zero);

    if(verbose_) 
    {
        forAll(diffWorkField.primitiveFieldRef(), cellI)
        {
            Utotal1 += diffWorkField.primitiveFieldRef()[cellI];
        }
    }

	diffusionTimeCount_[0] += particleCloud_.mesh().time().elapsedCpuTime() - t0;
    t0 = particleCloud_.mesh().time().elapsedCpuTime();
	
	while (diffusionRunTime_.loop())
    {
        if (diffusionRunTime_.timeIndex() == 1)
        {
            while (simple_.correctNonOrthogonal())
            {
                solve(
                    fvm::ddt(diffWorkField) 
                    - fvm::laplacian(DT, diffWorkField)
                );
            }
        }
        else
        {
            solve(
                    fvm::ddt(diffWorkField) 
                    - fvm::laplacian(DT, diffWorkField)
                );
        }
    }

    diffusionTimeCount_[1] += particleCloud_.mesh().time().elapsedCpuTime() - t0;
    t0 = particleCloud_.mesh().time().elapsedCpuTime();

    diffusionRunTime_.setTime(startTime,startIndex);

    vector Utotal2(vector::zero);

    if(verbose_) 
    {
        forAll(diffWorkField.primitiveFieldRef(), cellI)
        {
            Utotal2 += diffWorkField.primitiveFieldRef()[cellI]; //* alphas.primitiveFieldRef()[cellI];
        }
    }

	forAll(diffWorkField.primitiveFieldRef(), cellI)
    {
        if (alphas.primitiveFieldRef()[cellI] > ROOTVSMALL)
        {
            diffWorkField.primitiveFieldRef()[cellI] /=
                (alphas.primitiveFieldRef()[cellI]);
        }
    }
    // get data from working diffWorkField - will copy only values at new time
    fieldSrc.primitiveFieldRef() = diffWorkField.primitiveFieldRef();
    
    diffusionTimeCount_[0] += particleCloud_.mesh().time().elapsedCpuTime() - t0;
	
    fieldSrc.correctBoundaryConditions(); 

    if(verbose_)
    {   
        reduce(Utotal1,sumOp<vector>());
        reduce(Utotal2,sumOp<vector>());
        Info << "total Us*alphas before smoothing: " << Utotal1 << endl;
        Info << "total Us*alphas after smoothing: " << Utotal2 << endl;
        Info << "min/max(Us): " << min(fieldSrc) << tab << max(fieldSrc) << endl;
    }
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::diffusionSmoothing::smoothenReferenceField(volVectorField& fieldSrc) const   // for pure virtual function
{
	// do nothing
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
