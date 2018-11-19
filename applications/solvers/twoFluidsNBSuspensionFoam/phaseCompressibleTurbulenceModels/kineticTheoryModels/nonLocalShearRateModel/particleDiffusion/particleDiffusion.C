/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Author
    Federico Municchi, Purdue University (2018). All rights reserved.
\*---------------------------------------------------------------------------*/

#include "particleDiffusion.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
#include "fvc.H"
#include "twoPhaseSystem.H"
#include "slipFvPatchFields.H"
#include "kineticTheoryModel.H"
#include "phaseCompressibleTurbulenceModelFwd.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace nonLocalShearRateModels
{
    defineTypeNameAndDebug(particleDiffusion, 0);

    addToRunTimeSelectionTable
    (
        nonLocalShearRateModel,
        particleDiffusion,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::nonLocalShearRateModels::particleDiffusion::
particleDiffusion
(
    const dictionary& dict
)
:
    nonLocalShearRateModel(dict),
    coeffDict_(dict.subDict(typeName + "Coeffs")),
//    alpha_(readScalar(coeffDict_.lookup("alpha"))),
    Mus_(readScalar(coeffDict_.lookup("Mus"))),
    A_(readScalar(coeffDict_.lookup("A"))),
//    b_(readScalar(coeffDict_.lookup("b"))),
    nSteps_(readScalar(coeffDict_.lookup("nSteps")))
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::nonLocalShearRateModels::particleDiffusion::
~particleDiffusion()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::nonLocalShearRateModels::particleDiffusion::
gammaDotEff
(
   const phaseModel& phase,
   const volScalarField& Theta
) const
{
    tmp<volScalarField> tSR
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("tSR", phase.name()),
                Theta.time().timeName(),
                Theta.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            Theta.mesh(),
            dimensionedScalar("zero", dimless/dimTime, 0.0)
        )
    );

    volScalarField& SR = tSR.ref();

    //- First eveluate local shear rate
    const volVectorField& U = phase.U();
    volSymmTensorField shear = dev(symm(fvc::grad(U)));
    SR = sqrt(scalar(2)*(shear&&shear));

 /*   //- Then compute the yield coefficients based on the particle stresses
    volSymmTensorField R = -phase.turbulence().R();

    dimensionedScalar eps("eps",R.dimensions(),small);

    volScalarField P = -scalar(1./3.)*tr(R) + eps;

    R = R + P*symmTensor::I; 

    volScalarField Mu = sqrt(0.5*(R&&R))/P + small;

    //- Scale the shear rate with the yield coefficients
   // SR /= Mu;

    dimensionedScalar one("one",dimless,scalar(1));
    dimensionedScalar MUSdim("MUdim",Mu.dimensions(),Mus_);
       

    //- Cooperativity length
    tmp<volScalarField> zeta
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("Gamma", phase.name()),
                Theta.time().timeName(),
                Theta.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),

            A_
          * pow( 
                    (
                       one + Foam::pos(MUSdim-Mu))/(mag(Mu - MUSdim)),
                       alpha_
                    )
          * phase.d()    
        )
    );

    //- Particle mass
    volScalarField mp = 
        sqr(phase.d())*Foam::constant::mathematical::pi/scalar(6);

    //- Evaluate gloc
    volScalarField gloc = 
        Foam::pos(Mu-MUSdim)*(Mu-MUSdim)/(b_*Mu)*sqrt(mag(P)/mp);
*/
     dimensionedScalar MUSdim("MUdim",SR.dimensions(),Mus_);
    //- Diffusion coefficient
    volScalarField Gamma = SR*sqr(phase.d());    

    //- Solve diffusion equation explicitly (multi-step)

    dimensionedScalar subDeltaT = U.mesh().time().deltaT()/nSteps_; 

    for(label step=0;step<nSteps_;step++)
    {
        SR = SR + subDeltaT*fvc::laplacian(Gamma,SR);
        Gamma = SR*sqr(phase.d()); 
    }
    //- Re-scale shear rate
   // SR *= Mu;

    return tSR;

}



bool Foam::kineticTheoryModels::nonLocalShearRateModels::particleDiffusion::
read()
{
    coeffDict_ <<= dict_.subDict(typeName + "Coeffs");


    return true;
}


// ************************************************************************* //
