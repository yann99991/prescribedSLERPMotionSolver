# **Prescribed SLERP Motion Solver**

This folder contains the source code for the **OpenFOAM-ESI** version and two test cases to assist users in setting up simulations using the `prescribedSLERPMotionSolver`.

## **Overview**

The `prescribedSLERPMotionSolver` is a motion solver for OpenFOAM that allows for the definition of complex mesh motions using Spherical Linear Interpolation (SLERP). A detailed description of this code and its applications has been submitted to the OpenFOAM Journal for publication.

## **Test Cases**

Two motion types are included:

* **Harmonic** motion
* **Non-harmonic** motion

The geometry is that of a cylinder. The original 2D mesh was taken from the OpenFOAM tutorial case `oscillatingCylinder`, available via [this link](https://figshare.com/articles/presentation/OpenFOAM_advanced_training_Moving_meshes_rigid_body_motion_adaptive_mesh_refinement_and_overset_meshes/19310492).

The mesh has been **extruded in the *z*-direction** to produce a fully 3D domain, allowing the application of non-uniform deflections. A **U-bent** type deformation is applied in both cases, with the **mid-span section of the cylinder experiencing the maximum displacement**, while the edges of the cylinder remain undeformed. The deformation varies either harmonically or non-harmonically depending on the case.

---

### **Files Required**

To correctly set up and run each case, ensure the following files are included:

* `dynamicMeshDict` - defines the mesh motion solver and parameters
* Motion data files:

  * `motionData.dat` for harmonic motion
  * A series of `motionData_`*`TimeStep`*`.dat` files for non-harmonic motion
* `timeArray.dat` - required for the non-harmonic case to define motion timing
* `pointDisplacement` file in the `0/` directory - to specify displacement behaviour

Note: Since the `pointDisplacement` field is calculated by the solver for each mesh point, all boundaries of the domain can be assigned the type `calculated`. Additionally, for the non-harmonic motion case, the file `motionData_0.000000000.dat` must always be included, with all translational and rotational degrees of freedom set to zero.

---

### **Execution**

Each test case folder contains two `Allrun` scripts:

1. **`Allrun_moveDynamicMesh`** – Runs:

   * `checkMesh`
   * `moveDynamicMesh`:
     This is useful to verify mesh quality and observe the prescribed motion without solving the flow. The time-step size used here is larger, as fluid equations are not solved.
   * The `moveDynamicMesh` solver runs in **serial** and quickly generates the deformed mesh output at 0.1 s intervals.  

2. **`Allrun_pimpleFoam`** – Runs:

   * `checkMesh`
   * `pimpleFoam`:
     This solves the full case.
   * Default configuration runs in **parallel on 30 cores**, but this can be adjusted to suit your system and requirements.


---

### **Note**

These test cases are intended solely for educational or demonstration purposes, to illustrate the use of the `prescribedSLERPMotionSolver` in OpenFOAM. Further refinement may be necessary for use in physical or production-level simulations.


