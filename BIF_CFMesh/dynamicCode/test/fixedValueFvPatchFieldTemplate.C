/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "fixedValueFvPatchFieldTemplate.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "unitConversion.H"
//{{{ begin codeInclude

//}}} end codeInclude


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    // dynamicCode:
    // SHA1 = 3b0b44da1ee86de0ee54187b4bfd5a47bf4f9cb9
    //
    // unique function name that can be checked if the correct library version
    // has been loaded
    void test_3b0b44da1ee86de0ee54187b4bfd5a47bf4f9cb9(bool load)
    {
        if (load)
        {
            // code that can be explicitly executed after loading
        }
        else
        {
            // code that can be explicitly executed before unloading
        }
    }
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

makeRemovablePatchTypeField
(
    fvPatchVectorField,
    testFixedValueFvPatchVectorField
);


const char* const testFixedValueFvPatchVectorField::SHA1sum =
    "3b0b44da1ee86de0ee54187b4bfd5a47bf4f9cb9";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

testFixedValueFvPatchVectorField::
testFixedValueFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF)
{
    if (false)
    {
        Info<<"construct test sha1: 3b0b44da1ee86de0ee54187b4bfd5a47bf4f9cb9"
            " from patch/DimensionedField\n";
    }
}


testFixedValueFvPatchVectorField::
testFixedValueFvPatchVectorField
(
    const testFixedValueFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper)
{
    if (false)
    {
        Info<<"construct test sha1: 3b0b44da1ee86de0ee54187b4bfd5a47bf4f9cb9"
            " from patch/DimensionedField/mapper\n";
    }
}


testFixedValueFvPatchVectorField::
testFixedValueFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict)
{
    if (false)
    {
        Info<<"construct test sha1: 3b0b44da1ee86de0ee54187b4bfd5a47bf4f9cb9"
            " from patch/dictionary\n";
    }
}


testFixedValueFvPatchVectorField::
testFixedValueFvPatchVectorField
(
    const testFixedValueFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf)
{
    if (false)
    {
        Info<<"construct test sha1: 3b0b44da1ee86de0ee54187b4bfd5a47bf4f9cb9"
            " as copy\n";
    }
}


testFixedValueFvPatchVectorField::
testFixedValueFvPatchVectorField
(
    const testFixedValueFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF)
{
    if (false)
    {
        Info<<"construct test sha1: 3b0b44da1ee86de0ee54187b4bfd5a47bf4f9cb9 "
            "as copy/DimensionedField\n";
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

testFixedValueFvPatchVectorField::
~testFixedValueFvPatchVectorField()
{
    if (false)
    {
        Info<<"destroy test sha1: 3b0b44da1ee86de0ee54187b4bfd5a47bf4f9cb9\n";
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void testFixedValueFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (false)
    {
        Info<<"updateCoeffs test sha1: 3b0b44da1ee86de0ee54187b4bfd5a47bf4f9cb9\n";
    }

//{{{ begin code
    #line 36 "/hpc_atog/bsha219/OpenFoam-Tutorials/BIF_CFMesh/0/U.boundaryField.inlet"
const fvPatch& boundaryPatch = this->patch(); 
        const vectorField& Cf = boundaryPatch.Cf(); 

        vectorField& v = *this;
        scalar Uc = 0.05; // centerline velocity, U_max
	scalar R = 0.0084; // inlet radius

	forAll(Cf,faceI)
	{
	    scalar x = Cf[faceI].x(), z = Cf[faceI].z(), rSq = (x+0.1486715)*(x+0.1486715) + (z+0.203453)*(z+0.203453);
            if ( rSq < R )
             {
                v[faceI] = vector(0.2255*Uc*(1 - rSq/R/R), 0.9742*Uc*(1 - rSq/R/R), 0);
             };
            if ( rSq > R )
             {
               v[faceI] = vector(0, 0, 0);
             };
	}
//}}} end code

    this->fixedValueFvPatchField<vector>::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

