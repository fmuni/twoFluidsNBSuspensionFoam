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
    Foam::dragModels::MorrisBoulay

Description
    J.F. Morris & F. Boulay, 1999 Microstructure of strongly sheared 
    suspensions and its impact on rheology and diffusion.
    J. Rheol. 43, 1213–1237

SourceFiles
    MorrisBoulay.C
    
Author
    Federico Municchi, Purdue University (2018). All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef MorrisBoulay_H
#define MorrisBoulay_H

#include "dragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace dragModels
{

/*---------------------------------------------------------------------------*\
                           Class MorrisBoulay Declaration
\*---------------------------------------------------------------------------*/

class MorrisBoulay
:
    public dragModel
{
    
    //- Fitting parameter as in Donstov et al. (2013)
    dimensionedScalar alphaFit_;
    
public:

    //- Runtime type information
    TypeName("MorrisBoulay");


    // Constructors

        //- Construct from a dictionary and a phase pair
        MorrisBoulay
        (
            const dictionary& dict,
            const phasePair& pair,
            const bool registerObject
        );


    //- Destructor
    virtual ~MorrisBoulay();


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
