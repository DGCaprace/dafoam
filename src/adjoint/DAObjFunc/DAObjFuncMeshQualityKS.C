/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v3

\*---------------------------------------------------------------------------*/

#include "DAObjFuncMeshQualityKS.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(DAObjFuncMeshQualityKS, 0);
addToRunTimeSelectionTable(DAObjFunc, DAObjFuncMeshQualityKS, dictionary);
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

DAObjFuncMeshQualityKS::DAObjFuncMeshQualityKS(
    const fvMesh& mesh,
    const DAOption& daOption,
    const DAModel& daModel,
    const DAIndex& daIndex,
    const DAResidual& daResidual,
    const word objFuncName,
    const word objFuncPart,
    const dictionary& objFuncDict)
    : DAObjFunc(
        mesh,
        daOption,
        daModel,
        daIndex,
        daResidual,
        objFuncName,
        objFuncPart,
        objFuncDict)
{
    // Assign type, this is common for all objectives
    objFuncDict_.readEntry<word>("type", objFuncType_);

    objFuncConInfo_ = {};

    objFuncDict_.readEntry<scalar>("scale", scale_);

    objFuncDict_.readEntry<scalar>("coeffKS", coeffKS_);

    objFuncDict_.readEntry<word>("metric", metric_);

    // the polyMeshTools funcs are not fully AD in parallel, so the mesh quality
    // computed at the processor faces will not be properly back-propagate in AD
    // we can ignore the mesh quality at proc patches if includeProcPatches = False (default)
    includeProcPatches_ = objFuncDict_.lookupOrDefault<label>("includeProcPatches", 0);
    includeFaceList_.setSize(mesh_.nFaces(), 1);
    if (!includeProcPatches_)
    {
        label nInterF = mesh_.nInternalFaces();
        label faceCounter = 0;
        forAll(mesh.boundaryMesh(), patchI)
        {
            if (mesh.boundaryMesh()[patchI].type() == "processor")
            {
                forAll(mesh.boundaryMesh()[patchI], faceI)
                {
                    includeFaceList_[nInterF + faceCounter] = 0;
                    faceCounter += 1;
                }
            }
            else
            {
                faceCounter += mesh.boundaryMesh()[patchI].size();
            }
        }
    }

    if (daOption.getOption<label>("debug"))
    {
        Info << "includeFaceList " << includeFaceList_ << endl;
    }

}

/// calculate the value of objective function
void DAObjFuncMeshQualityKS::calcObjFunc(
    const labelList& objFuncFaceSources,
    const labelList& objFuncCellSources,
    scalarList& objFuncFaceValues,
    scalarList& objFuncCellValues,
    scalar& objFuncValue)
{
    /*
    Description:
        Calculate the mesh quality and aggregate with the KS function
        e.g., if metric is the faceSkewness, the objFunc value will be the
        approximated max skewness

    Input:
        objFuncFaceSources: List of face source (index) for this objective
    
        objFuncCellSources: List of cell source (index) for this objective

    Output:
        objFuncFaceValues: the discrete value of objective for each face source (index). 
        This  will be used for computing df/dw in the adjoint.
    
        objFuncCellValues: the discrete value of objective on each cell source (index). 
        This will be used for computing df/dw in the adjoint.
    
        objFuncValue: the sum of objective, reduced across all processors and scaled by "scale"
    */

    // initialize objFunValue
    objFuncValue = 0.0;

    if (metric_ == "faceOrthogonality")
    {
        // faceOrthogonality ranges from 0 to 1 and the nonOrthoAngle is acos(faceOrthogonality)
        const scalarField faceOrthogonality(
            polyMeshTools::faceOrthogonality(
                mesh_,
                mesh_.faceAreas(),
                mesh_.cellCentres()));

        // calculate the KS mesh quality
        forAll(faceOrthogonality, faceI)
        {
            if (includeFaceList_[faceI] == 1)
            {
                objFuncValue += exp(coeffKS_ * faceOrthogonality[faceI]);
            }

            if (objFuncValue > 1e200)
            {
                FatalErrorIn(" ") << "KS function summation term too large! "
                                  << "Reduce coeffKS! " << abort(FatalError);
            }
        }
    }
    else if (metric_ == "nonOrthoAngle")
    {
        const scalarField faceOrthogonality(
            polyMeshTools::faceOrthogonality(
                mesh_,
                mesh_.faceAreas(),
                mesh_.cellCentres()));

        // Face based non ortho angle
        scalarField nonOrthoAngle = faceOrthogonality;
        forAll(faceOrthogonality, faceI)
        {
            scalar val = faceOrthogonality[faceI];
            // bound it to less than 1.0 - 1e-6. We can't let val = 1
            // because its derivative will be divided by zero
            scalar boundV = 1.0 - 1e-6;
            if (val > boundV)
            {
                val = boundV;
            }
            if (val < -boundV)
            {
                val = -boundV;
            }
            // compute non ortho angle
            scalar angleRad = acos(val);
            // convert rad to degree
            scalar pi = constant::mathematical::pi;
            scalar angleDeg = angleRad * 180.0 / pi;
            nonOrthoAngle[faceI] = angleDeg;
        }

        // calculate the KS mesh quality
        forAll(nonOrthoAngle, faceI)
        {
            if (includeFaceList_[faceI] == 1)
            {
                objFuncValue += exp(coeffKS_ * nonOrthoAngle[faceI]);
            }

            if (objFuncValue > 1e200)
            {
                FatalErrorIn(" ") << "KS function summation term too large! "
                                  << "Reduce coeffKS! " << abort(FatalError);
            }
        }
    }
    else if (metric_ == "faceSkewness")
    {
        const scalarField faceSkewness(
            polyMeshTools::faceSkewness(
                mesh_,
                mesh_.points(),
                mesh_.faceCentres(),
                mesh_.faceAreas(),
                mesh_.cellCentres()));

        // calculate the KS mesh quality
        forAll(faceSkewness, faceI)
        {
            if (includeFaceList_[faceI] == 1)
            {
                objFuncValue += exp(coeffKS_ * faceSkewness[faceI]);
            }

            if (objFuncValue > 1e200)
            {
                FatalErrorIn(" ") << "KS function summation term too large! "
                                  << "Reduce coeffKS! " << abort(FatalError);
            }
        }
    }
    else
    {
        FatalErrorIn(" ") << "metric not valid! "
                          << "Options: faceOrthogonality, nonOrthoAngle, or faceSkewness "
                          << abort(FatalError);
    }

    // need to reduce the sum of force across all processors
    reduce(objFuncValue, sumOp<scalar>());

    objFuncValue = log(objFuncValue) / coeffKS_;

    return;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
