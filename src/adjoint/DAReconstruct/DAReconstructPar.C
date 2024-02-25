/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v3

\*---------------------------------------------------------------------------*/


#include "DAReconstructPar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Constructors
DAReconstructPar::DAReconstructPar(
    const DAOption& daOption)
    : daOption_(daOption)
{

    // fvMesh& meshNew = const_cast<fvMesh&>(mesh);


    // maxNonOrth_ = daOption_.getSubDictOption<scalar>("checkMeshThreshold", "maxNonOrth");
    // maxIncorrectlyOrientedFaces_ =
    //     daOption_.getSubDictOption<label>("checkMeshThreshold", "maxIncorrectlyOrientedFaces");
    
 
    Info << "Hellow from the DAReconstructPar constructor: " << endl;
    // Info << "maxNonOrth: " << maxNonOrth_ << endl;
    
    // word surfaceFormat = "vtk";
    // surfWriter.reset(surfaceWriter::New(surfaceFormat));
    // setWriter.reset(writer<scalar>::New(vtkSetWriter<scalar>::typeName));
}

DAReconstructPar::~DAReconstructPar()
{
}

label DAReconstructPar::run() const
{
    /*
    Description:
        Run checkMesh and return meshOK
    
    Output:
        meshOK: 1 means quality passes
    */

#if !defined(CODI_AD_FORWARD) & !defined(CODI_AD_REVERSE)    

    label reconstructOK = 1;

    Info << "reconstructPar" << endl;
    // Info << "reconstructPar for time = " << runTime.timeName() << endl;

    /*
    //We always reconstruct all regions
    wordList regionNames = regionProperties(runTime).names();
    wordList regionDirs = regionNames;

    Info<< "Reconstructing all regions in regionProperties" << nl
        << "    " << flatOutput(regionNames) << nl << endl;

    // Determine the processor count
    label nProcs = fileHandler().nProcs(args.path(), regionDirs[0]);

    // Warn fileHandler of number of processors
    const_cast<fileOperation&>(fileHandler()).setNProcs(nProcs);

    // Create the processor databases
    PtrList<Time> databases(nProcs);

    forAll(databases, proci)
    {
        databases.set
        (
            proci,
            new Time
            (
                Time::controlDictName,
                args.rootPath(),
                args.caseName()/("processor" + Foam::name(proci))
            )
        );
    }

    // TODO: tHE time should be the current time

    // Set all times on processor meshes equal to reconstructed mesh
    forAll(databases, proci)
    {
        databases[proci].setTime(runTime);
    }

    //For all regions:
    forAll(regionNames, regioni)
    {
        const word& regionName = regionNames[regioni];
        const word& regionDir = regionDirs[regioni];


        //DG: NO NEED TO READ THE MESH SINCE WE HAVE IT IN MEM
        // fvMesh mesh
        // (
        //     IOobject
        //     (
        //         regionName,
        //         runTime.timeName(),
        //         runTime,
        //         Foam::IOobject::MUST_READ
        //     )
        // );


        // Read all meshes and addressing to reconstructed mesh
        processorMeshes procMeshes(databases, regionName);


        //DG: THE MESH IS UP TO DATE BY CONSTRUCTION. SKIPPING THIS
        // // Check if any new meshes need to be read.
        // fvMesh::readUpdateState meshStat = mesh.readUpdate();
        // fvMesh::readUpdateState procStat = procMeshes.readUpdate();
        //
        // if (procStat == fvMesh::POINTS_MOVED)
        // {
        //     // Reconstruct the points for moving mesh cases and write
        //     // them out
        //     procMeshes.reconstructPoints(mesh);
        // }
        // else if (meshStat != procStat)
        // {
        //     WarningInFunction
        //         << "readUpdate for the reconstructed mesh:"
        //         << meshStat << nl
        //         << "readUpdate for the processor meshes  :"
        //         << procStat << nl
        //         << "These should be equal or your addressing"
        //         << " might be incorrect."
        //         << " Please check your time directories for any "
        //         << "mesh directories." << endl;
        // }

        // Get list of objects from processor0 database
        IOobjectList objects
        (
            procMeshes.meshes()[0],
            databases[0].timeName()
        );


        //dofields: always true here
        {
            // If there are any FV fields, reconstruct them
            Info<< "Reconstructing FV fields" << nl << endl;

            fvFieldReconstructor reconstructor
            (
                mesh,
                procMeshes.meshes(),
                procMeshes.faceProcAddressing(),
                procMeshes.cellProcAddressing(),
                procMeshes.boundaryProcAddressing()
            );

            reconstructor.reconstructFvVolumeInternalFields<scalar>
            (
                objects,
                selectedFields
            );
            reconstructor.reconstructFvVolumeInternalFields<vector>
            (
                objects,
                selectedFields
            );
            reconstructor.reconstructFvVolumeInternalFields<sphericalTensor>
            (
                objects,
                selectedFields
            );
            reconstructor.reconstructFvVolumeInternalFields<symmTensor>
            (
                objects,
                selectedFields
            );
            reconstructor.reconstructFvVolumeInternalFields<tensor>
            (
                objects,
                selectedFields
            );

            reconstructor.reconstructFvVolumeFields<scalar>
            (
                objects,
                selectedFields
            );
            reconstructor.reconstructFvVolumeFields<vector>
            (
                objects,
                selectedFields
            );
            reconstructor.reconstructFvVolumeFields<sphericalTensor>
            (
                objects,
                selectedFields
            );
            reconstructor.reconstructFvVolumeFields<symmTensor>
            (
                objects,
                selectedFields
            );
            reconstructor.reconstructFvVolumeFields<tensor>
            (
                objects,
                selectedFields
            );

            reconstructor.reconstructFvSurfaceFields<scalar>
            (
                objects,
                selectedFields
            );
            reconstructor.reconstructFvSurfaceFields<vector>
            (
                objects,
                selectedFields
            );
            reconstructor.reconstructFvSurfaceFields<sphericalTensor>
            (
                objects,
                selectedFields
            );
            reconstructor.reconstructFvSurfaceFields<symmTensor>
            (
                objects,
                selectedFields
            );
            reconstructor.reconstructFvSurfaceFields<tensor>
            (
                objects,
                selectedFields
            );

            if (reconstructor.nReconstructed() == 0)
            {
                Info<< "No FV fields" << nl << endl;
            }
        }

        //dofields: always true here
        {
            Info<< "Reconstructing point fields" << nl << endl;

            const pointMesh& pMesh = pointMesh::New(mesh);
            PtrList<pointMesh> pMeshes(procMeshes.meshes().size());

            forAll(pMeshes, proci)
            {
                pMeshes.set
                (
                    proci,
                    new pointMesh(procMeshes.meshes()[proci])
                );
            }

            pointFieldReconstructor reconstructor
            (
                pMesh,
                pMeshes,
                procMeshes.pointProcAddressing(),
                procMeshes.boundaryProcAddressing()
            );

            reconstructor.reconstructFields<scalar>
            (
                objects,
                selectedFields
            );
            reconstructor.reconstructFields<vector>
            (
                objects,
                selectedFields
            );
            reconstructor.reconstructFields<sphericalTensor>
            (
                objects,
                selectedFields
            );
            reconstructor.reconstructFields<symmTensor>
            (
                objects,
                selectedFields
            );
            reconstructor.reconstructFields<tensor>
            (
                objects,
                selectedFields
            );

            if (reconstructor.nReconstructed() == 0)
            {
                Info<< "No point fields" << nl << endl;
            }
        }


        // skipping Lagrangian fields
        // skipping FA fiels
        // skipping sets
        // skipping refinement data


        // If there is a "uniform" directory in the time region
        // directory copy from the master processor
        {
            fileName uniformDir0
            (
                fileHandler().filePath
                (
                    databases[0].timePath()/regionDir/"uniform"
                )
            );

            if (!uniformDir0.empty() && fileHandler().isDir(uniformDir0))
            {
                fileHandler().cp(uniformDir0, runTime.timePath()/regionDir);
            }
        }

        // For the first region of a multi-region case additionally
        // copy the "uniform" directory in the time directory
        if (regioni == 0 && regionDir != word::null)
        {
            fileName uniformDir0
            (
                fileHandler().filePath
                (
                    databases[0].timePath()/"uniform"
                )
            );

            if (!uniformDir0.empty() && fileHandler().isDir(uniformDir0))
            {
                fileHandler().cp(uniformDir0, runTime.timePath());
            }
        }

        // if (nFailedChecks)
        // {
        //     Info << "\nFailed " << nFailedChecks << " mesh checks.\n"
        //         << endl;
        //     meshOK = 0;
        // }
        // else
        // {
        //     Info << "\nMesh OK.\n"
        //         << endl;
        // }
    }

    */

    return reconstructOK;
#else
    return 1;
#endif    
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //