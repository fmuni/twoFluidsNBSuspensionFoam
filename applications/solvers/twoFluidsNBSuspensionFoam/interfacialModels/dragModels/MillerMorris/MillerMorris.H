/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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
    Foam::dragModels::MillerMorris

Description
    Miller M., Morris F., 2006, Normal stress-driven migration and axial 
    development in pressure-driven flow of concentrated suspensions.
    Journal of Non-Newtonian Fluid Mechanics 135, 149-165

SourceFiles
    MillerMorris.C
    
Author
    Federico Municchi, Purdue University (2018). All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef MillerMorris_H
#define MillerMorris_H

#include "dragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace dragModels
{

/*---------------------------------------------------------------------------*\
                           Class MillerMorris Declaration
\*---------------------------------------------------------------------------*/

class MillerMorris
:
    public dragModel
{
    
    //- Fitting parameter as in Donstov et al. (2013)
    dimensionedScalar alphaFit_;
    
public:

    //- Runtime type information
    TypeName("MillerMorris");


    // Constructors

        //- Construct from a dictionary and a phase pair
        MillerMorris
        (
            const dictionary& dict,
            const phasePair& pair,
            const bool registerObject
        );


    //- Destructor
    virtual ~MillerMorris();


    // Member Functions

        //- Drag coefficient
        virtual tmp<volScalarField> CdRe() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dragModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
