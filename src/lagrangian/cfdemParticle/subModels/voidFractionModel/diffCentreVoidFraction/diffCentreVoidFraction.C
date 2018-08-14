/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright 2009-2012 JKU Linz
                                Copyright 2012-     DCS Computing GmbH, Linz
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

Author 
    Linhan Ge, University of Newcastle
    linhan.ge@gmail.com
\*---------------------------------------------------------------------------*/

#include "error.H"

#include "diffCentreVoidFraction.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(diffCentreVoidFraction, 0);

addToRunTimeSelectionTable
(
    voidFractionModel,
    diffCentreVoidFraction,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
diffCentreVoidFraction::diffCentreVoidFraction
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    voidFractionModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props"))
{
    checkWeightNporosity(propsDict_);
    if(porosity()!=1) FatalError << "porosity not used in diffCentreVoidFraction" << abort(FatalError);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

diffCentreVoidFraction::~diffCentreVoidFraction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void diffCentreVoidFraction::setvoidFraction(double** const& mask,double**& voidfractions,double**& particleWeights,double**& particleVolumes,double**& particleV) const
{
    reAllocArrays();

    scalar radius(-1);
    scalar volume(0);
    scalar cellVol(0);
    scalar scaleVol= weight();

    for(int index=0; index< particleCloud_.numberOfParticles(); index++)
    {
            particleWeights[index][0]=0;
            cellsPerParticle_[index][0]=1;

            label cellI = particleCloud_.cellIDs()[index][0];

            if (cellI >= 0)  // particel centre is in domain
            {
                cellVol = particlefractionNext_.mesh().V()[cellI];    // voidfractionNext_ is a volScalarField, the member function mesh is??
                radius = particleCloud_.radius(index);
                volume = 4.188790205*radius*radius*radius*scaleVol;

                // store volume for each particle
                particleVolumes[index][0] = volume/cellVol;           // in diffusion method, this should be volume fraction,pass to setVoidFraction
                particleV[index][0] = volume;

                particlefractionNext_[cellI] += volume/cellVol;				

                // store cellweight for each particle  - this should not live here
                particleWeights[index][0] = 1;
            }
    }
    particlefractionNext_.correctBoundaryConditions();

    // bring voidfraction from Eulerian Field to particle array
    for(int index=0; index< particleCloud_.numberOfParticles(); index++)
    {
        // label cellID = particleCloud_.cellIDs()[index][0];

        voidfractions[index][0] = -1;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
