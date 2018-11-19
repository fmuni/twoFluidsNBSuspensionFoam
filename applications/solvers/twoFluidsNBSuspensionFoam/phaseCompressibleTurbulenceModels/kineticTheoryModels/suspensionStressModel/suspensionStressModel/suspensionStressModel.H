/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::kineticTheoryModels::suspensionStressModel

SourceFiles
    suspensionStressModel.C

\*---------------------------------------------------------------------------*/

#ifndef suspensionStressModel_H
#define suspensionStressModel_H

#include "dictionary.H"
#include "volFields.H"
#include "dimensionedTypes.H"
#include "runTimeSelectionTables.H"
#include "phaseModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{

/*---------------------------------------------------------------------------*\
                           Class suspensionStressModel Declaration
\*---------------------------------------------------------------------------*/

class suspensionStressModel
{
    // Private member functions

        //- Disallow default bitwise copy construct
        suspensionStressModel(const suspensionStressModel&);

        //- Disallow default bitwise assignment
        void operator=(const suspensionStressModel&);


protected:

    // Protected data

        //- Reference to higher-level dictionary for re-read
        const dictionary& dict_;


public:

    //- Runtime type information
    TypeName("suspensionStressModel");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        suspensionStressModel,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );


    // Constructors

        //- Construct from components
        suspensionStressModel(const dictionary& dict);


    // Selectors

        static autoPtr<suspensionStressModel> New
        (
            const dictionary& dict
        );


    //- Destructor
    virtual ~suspensionStressModel();


    // Member Functions

        // Return directional anisotropic tensor
        virtual tmp<volSymmTensorField> Q
        (
            const volVectorField& U,
            const vector& lambda
        ) const;

        virtual tmp<volSymmTensorField> Q
        (
            const volVectorField& U,
            const volScalarField& lambda1,
            const volScalarField& lambda2,
            const volScalarField& lambda3
        ) const;

        virtual tmp<volSymmTensorField> computeQ
        (
            const phaseModel& phase
        ) const = 0;

        virtual tmp<volVectorField> shearForce
        (
            const phaseModel& phase,
            const dimensionedScalar& alphaMax,
            const volScalarField& gammaDotEff,
            const volSymmTensorField& fabT
        ) const = 0;

        virtual tmp<volSymmTensorField> DSigmaDalpha
        (
            const phaseModel& phase,
            const dimensionedScalar& alphaMax,
            const volScalarField& gammaDotEff,
            const volSymmTensorField& fabT
        ) const = 0;

        virtual tmp<volSymmTensorField> SigmaS
        (
            const phaseModel& phase,
            const dimensionedScalar& alphaMax,
            const volScalarField& gammaDotEff,
            const volSymmTensorField& fabT
        ) const = 0;

        virtual tmp<volScalarField> nu
        (
            const phaseModel& phase,
            const dimensionedScalar& alphaMax
        ) const = 0;

        virtual bool read() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace kineticTheoryModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
