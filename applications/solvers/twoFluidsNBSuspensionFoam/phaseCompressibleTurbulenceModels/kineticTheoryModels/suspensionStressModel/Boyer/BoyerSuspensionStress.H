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

Class
    Foam::kineticTheoryModels::suspensionStressModels::Boyer

Description
    F. Boyer, E. Guazzelli & O. Poliquen, 2011, Unifying suspension
    and granular rheology. Phys. Rev. Lett. 107, 188301.

SourceFiles
    BoyerSuspensionStress.C

Author
    Federico Municchi, Purdue University (2018). All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef BoyerSuspensionStress_H
#define BoyerSuspensionStress_H

#include "suspensionStressModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace suspensionStressModels
{

/*---------------------------------------------------------------------------*\
               Class Boyer Declaration
\*---------------------------------------------------------------------------*/

class Boyer
:
    public suspensionStressModel
{
    // Private data

        dictionary coeffDict_;

        //- Material constant for suspension normal stress
        dimensionedScalar mu1_;

        //- Material constant for suspension normal stress
        dimensionedScalar mu2_;

        //- Material constant for suspension normal stress
        dimensionedScalar I0_;

        //- Anisotropy coefficients for deviatoric pressure
        dimensionedVector lambda_;



public:

    //- Runtime type information
    TypeName("Boyer");


    // Constructors

        //- Construct from components
        Boyer(const dictionary& dict);


    //- Destructor
    virtual ~Boyer();


    // Member functions
        
        virtual tmp<volSymmTensorField> computeQ
        (
            const phaseModel& phase
        ) const ;

        virtual tmp<volVectorField> shearForce
        (
            const phaseModel& phase,
            const dimensionedScalar& alphaMax,
            const volScalarField& gammaDotEff,
            const volSymmTensorField& fabT
        ) const ;

        virtual tmp<volSymmTensorField> DSigmaDalpha
        (
            const phaseModel& phase,
            const dimensionedScalar& alphaMax,
            const volScalarField& gammaDotEff,
            const volSymmTensorField& fabT
        )  const;

        virtual tmp<volSymmTensorField> SigmaS
        (
            const phaseModel& phase,
            const dimensionedScalar& alphaMax,
            const volScalarField& gammaDotEff,
            const volSymmTensorField& fabT
        ) const;        

        virtual tmp<volScalarField> nu
        (
            const phaseModel& phase,
            const dimensionedScalar& alphaMax
        ) const ;

        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace suspensionStressModels
} // End namespace kineticTheoryModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
