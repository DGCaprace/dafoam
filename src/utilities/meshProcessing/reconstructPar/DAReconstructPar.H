/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v3

    Description:
        Check mesh quality, modified from the checkMesh utility

\*---------------------------------------------------------------------------*/

#ifndef DAReconstructPar_H
#define DAReconstructPar_H

#include "DAOption.H"
#include "Time.H"
#include "fvMesh.H"
#include "IOdictionary.H"

// #include "processorMeshes.H" // <---------------------------?
// #include "regionProperties.H" // <---------------------------?
// #include "fvFieldReconstructor.H" 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class DAReconstructPar Declaration
\*---------------------------------------------------------------------------*/

class DAReconstructPar
{

private:
    /// Disallow default bitwise copy construct
    DAReconstructPar(const DAReconstructPar&);

    /// Disallow default bitwise assignment
    void operator=(const DAReconstructPar&);

    /// DAOption object
    const DAOption& daOption_;


    // DO I NEED THIS??
    // #include "setRootCase.H"
    // #include "createTime.H"

    /// Time object
    const Time& runTime;

    /// fvMesh object
    const fvMesh& mesh;

public:
    // Constructors

    /// Construct from components
    DAReconstructPar(
        const DAOption& daOption,
        const Time& runTime1,
        const fvMesh& mesh1);

    /// Destructor
    virtual ~DAReconstructPar();

    /// run reconstructPar and return reconstructOK
    label run() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
