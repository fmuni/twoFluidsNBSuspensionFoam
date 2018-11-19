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
    Foam::kineticTheoryModels::suspensionStressModels::DLL

Description
    Dbouk, T., Lobry, L., & Lemaire, E. (2013). Normal stresses in concentrated
    non-Brownian suspensions. Journal of Fluid Mechanics, 715, 239-272. 
    doi:10.1017/jfm.2012.516

SourceFiles
    DboukLobryLemaireSuspensionStress.C

Author
    Federico Municchi, Purdue University (2018). All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef DLLSuspensionStress_H
#define DLLSuspensionStress_H

#include "suspensionStressModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace suspensionStressModels
{

/*---------------------------------------------------------------------------*\
               Class DLL Declaration
\*---------------------------------------------------------------------------*/

class DLL
:
    public suspensionStressModel
{
    // Private data

        dictionary coeffDict_;

        //- Material constant for suspension normal stress
        dimensionedScalar a_;

        //- Material constant for suspension normal stress
        dimensionedScalar b_;

        //- Material constant for suspension normal stress
        dimensionedScalar c_;

        //- Material constant for suspension normal stress
        dimensionedScalar KN_;

        //- Use fixed lambdas
        Switch fixedLambdas_;




public:

    //- Runtime type information
    TypeName("DLL");


    // Constructors

        //- Construct from components
        DLL(const dictionary& dict);


    //- Destructor
    virtual ~DLL();


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
