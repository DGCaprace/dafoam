#!/usr/bin/env python
"""
Run Python tests for optimization integration
"""

from mpi4py import MPI
import os
import numpy as np
from testFuncs import *

import openmdao.api as om
from mphys.multipoint import Multipoint
from dafoam.mphys import DAFoamBuilder, OptFuncs
from mphys.scenario_aerodynamic import ScenarioAerodynamic
from pygeo.mphys import OM_DVGEOCOMP
from pygeo import geo_utils

gcomm = MPI.COMM_WORLD

os.chdir("./input/NACA0012")
if gcomm.rank == 0:
    os.system("rm -rf 0 processor*")
    os.system("cp -r 0.incompressible 0")
    os.system("cp -r system.incompressible system")
    os.system("cp -r constant/turbulenceProperties.safv3 constant/turbulenceProperties")

# aero setup
U0 = 10.0
p0 = 0.0
A0 = 0.1
aoa0 = 3.0
LRef = 1.0

daOptions = {
    "designSurfaces": ["wing"],
    "solverName": "DASimpleFoam",
    "primalMinResTol": 1.0e-10,
    "primalMinResTolDiff": 1e4,
    "primalBC": {
        "U0": {"variable": "U", "patches": ["inout"], "value": [U0, 0.0, 0.0]},
        "p0": {"variable": "p", "patches": ["inout"], "value": [p0]},
        "useWallFunction": True,
        "transport:nu": 1.5e-5,
    },
    "objFunc": {
        "CL": {
            "part1": {
                "type": "force",
                "source": "patchToFace",
                "patches": ["wing"],
                "directionMode": "normalToFlow",
                "alphaName": "aoa",
                "scale": 1.0 / (0.5 * U0 * U0 * A0),
                "addToAdjoint": True,
            }
        },
    },
    "adjEqnOption": {
        "gmresRelTol": 1.0e-8,
        "pcFillLevel": 1,
        "jacMatReOrdering": "rcm",
    },
    "normalizeStates": {"U": U0, "p": U0 * U0 / 2.0, "phi": 1.0, "nuTilda": 1e-3},
    "adjPartDerivFDStep": {"State": 1e-6},
    "designVar": {
        "twist": {"designVarType": "FFD"},
        "aoa": {"designVarType": "AOA", "patches": ["inout"], "flowAxis": "x", "normalAxis": "y"},
        "alphaPorosity": {"designVarType": "Field", "fieldName": "alphaPorosity", "fieldType": "scalar"},
    },
}

meshOptions = {
    "gridFile": os.getcwd(),
    "fileType": "OpenFOAM",
    # point and normal for the symmetry plane
    "symmetryPlanes": [[[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]], [[0.0, 0.0, 0.1], [0.0, 0.0, 1.0]]],
}


class Top(Multipoint):
    def setup(self):
        dafoam_builder = DAFoamBuilder(daOptions, meshOptions, scenario="aerodynamic")
        dafoam_builder.initialize(self.comm)

        ################################################################################
        # MPHY setup
        ################################################################################

        # ivc to keep the top level DVs
        self.add_subsystem("dvs", om.IndepVarComp(), promotes=["*"])

        # create the mesh and cruise scenario because we only have one analysis point
        self.add_subsystem("mesh", dafoam_builder.get_mesh_coordinate_subsystem())

        # add the geometry component, we dont need a builder because we do it here.
        self.add_subsystem("geometry", OM_DVGEOCOMP(file="FFD/wingFFD.xyz", type="ffd"))

        self.mphys_add_scenario("cruise", ScenarioAerodynamic(aero_builder=dafoam_builder))

        self.connect("mesh.x_aero0", "geometry.x_aero_in")
        self.connect("geometry.x_aero0", "cruise.x_aero")

    def configure(self):

        self.cruise.aero_post.mphys_add_funcs()

        # create geometric DV setup
        points = self.mesh.mphys_get_surface_mesh()

        # add pointset
        self.geometry.nom_add_discipline_coords("aero", points)

        # geometry setup

        # Create reference axis
        def aoa(val, DASolver):
            aoa = val[0] * np.pi / 180.0
            U = [float(U0 * np.cos(aoa)), float(U0 * np.sin(aoa)), 0]
            # we need to update the U value only
            DASolver.setOption("primalBC", {"U0": {"value": U}})
            DASolver.updateDAOption()

        # pass this aoa function to the cruise group
        self.cruise.coupling.solver.add_dv_func("aoa", aoa)
        self.cruise.aero_post.add_dv_func("aoa", aoa)

        def alphaPorosity(val, DASolver):
            for idxI, v in enumerate(val):
                DASolver.setFieldValue4LocalCellI(b"alphaPorosity", v, idxI)
                DASolver.updateBoundaryConditions(b"alphaPorosity", b"scalar")

        self.cruise.coupling.solver.add_dv_func("alphaPorosity", alphaPorosity)
        self.cruise.aero_post.add_dv_func("alphaPorosity", alphaPorosity)

        # Set up global design variables. We dont change the root twist
        nRefAxPts = self.geometry.nom_addRefAxis(name="wingAxis", xFraction=0.25, alignIndex="k")

        def twist(val, geo):
            for i in range(nRefAxPts):
                geo.rot_z["wingAxis"].coef[i] = -val[0]

        self.geometry.nom_addGlobalDV(dvName="twist", value=np.array([0]), func=twist)

        # add dvs to ivc and connect
        self.dvs.add_output("twist", val=np.array([0]))
        self.dvs.add_output("aoa", val=np.array([aoa0]))
        nLocalCells = self.cruise.coupling.DASolver.solver.getNLocalCells()
        self.dvs.add_output("alphaPorosity", val=np.zeros(nLocalCells, dtype="d"))

        self.connect("twist", "geometry.twist")
        self.connect("aoa", "cruise.aoa")
        self.connect("alphaPorosity", "cruise.alphaPorosity")

        # define the design variables
        self.add_design_var("twist", lower=-10.0, upper=10.0, scaler=1.0)
        self.add_design_var("aoa", lower=-10.0, upper=10.0, scaler=1.0)
        self.add_design_var("alphaPorosity", lower=0, upper=1e-4, scaler=1.0)

        # add constraints and the objective
        self.add_constraint("cruise.aero_post.CL", equals=0.3, scaler=1.0)


prob = om.Problem(reports=None)
prob.model = Top()

prob.driver = om.pyOptSparseDriver()
prob.driver.options["optimizer"] = "IPOPT"
prob.driver.opt_settings = {
    "tol": 1.0e-7,
    "max_iter": 50,
    "output_file": "opt_IPOPT.out",
    "constr_viol_tol": 1.0e-7,
    "mu_strategy": "adaptive",
    "limited_memory_max_history": 26,
    "nlp_scaling_method": "gradient-based",
    "alpha_for_y": "full",
    "recalc_y": "yes",
    "print_level": 5,
    "acceptable_tol": 1.0e-7,
}
prob.driver.options["debug_print"] = ["nl_cons", "objs", "desvars"]


prob.add_recorder(om.SqliteRecorder("cases.sql"))
prob.recording_options["includes"] = []
prob.recording_options["record_objectives"] = True
prob.recording_options["record_constraints"] = True

prob.setup(mode="rev")
om.n2(prob, show_browser=False, outfile="mphys_aero.html")

optFuncs = OptFuncs(daOptions, prob)

prob.run_model()
totals = prob.compute_totals()

alphaNorm = np.linalg.norm(totals[("cruise.aero_post.functionals.CL", "dvs.alphaPorosity")].flatten())
alphaNormSum = gcomm.allreduce(alphaNorm, op=MPI.SUM)

if gcomm.rank == 0:
    objFuncDict = {}
    objFuncDict["CL"] = prob.get_val("cruise.aero_post.CL")
    derivDict = {}
    derivDict["CL"] = {}
    derivDict["CL"]["alpha"] = [alphaNormSum]
    derivDict["CL"]["twist"] = totals[("cruise.aero_post.functionals.CL", "dvs.twist")][0]
    derivDict["CL"]["aoa"] = totals[("cruise.aero_post.functionals.CL", "dvs.aoa")][0]
    reg_write_dict(objFuncDict, 1e-6, 1e-8)
    reg_write_dict(derivDict, 1e-4, 1e-6)
