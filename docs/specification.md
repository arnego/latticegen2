# latticegen2 Specification

What we currently know about the latticegen2 project. Where we don't know yet we write [TODO: needs decision] rather
than leaving it blank, so gaps are visible instead of silently assumed.
When a feature or characteristic of this project has been proposed by claude, it must clearly state so by using [TODO: proposed] rather than leaving it blank, so the user can control what enters the specification.
Do not implement something that is tagged with [TODO: needs decision] or [TODO: proposed]. 

---

## 1. Purpose & Scope

**Goal:**
The script must generate and output a parameterized lattice geometry based on user input that fits exactly within the users boundry geometry. 
The script must use a highly optimized and parameterized generation algorithm that can taking the running hardware into account in order to ensure minimum duration to output, and good stability of the run-time system.

**Primary output:** 
A single watertight STEP representing a lattice core, filling volume defined by the solid body of the input STEP geometry, with boundry against the surfaces of the input STEP geometry, placed within the same coordinate system as the input STEP geometry. The STEP file may contain multiple bodies if the input geometry cuts rods off in such a way that some rods become floating islands disconnected from the rest.

**Secondary output:** 
Run data from the script including:
  - The runs input parameters 
  - Date and time of run start
  - Duration from start to completion
  - Run characteristics (number of tiles, number of parallel threads per stage of the generation procedure, etc)
  - Maximum memory usage

---

## 2. Deployment Target & Constraints

- **Runtime environment:** Windows 11 offline workstation and Linux command line
- **Julia version:** 1.10
- **Offline requirement:** Package must run with **zero network access**. List anything that currently requires network (package managers, license checks, telemetry) so it can be eliminated or vendored.
- **Packaging form:** TODO
- **Target machine specs / limits:** Main development system: 32 GB RAM, 6 cure CPU, Nvidia RTX 3080 GPU , disk space for intermediate meshes
RAM and CPU cores should optinally be provided as input parameters, and optimization parameters should be determined automatically thereafter. If RAM and CPU cores are not provided, the optimization parameters must be provided explicitly instead. Script should run at a priority of one below normal to usability of the workstation.
- **Allowed third-party libraries:** Must be compatible with the target OS/arch. License text must be obtained and put into /licenses folder, and @/licenses/libraries.md must be updated with the cross reference between the library used and the corresponding license text file valid for that library.
- **License constraints:** TBD

---

## 3. Command-Line Interface

Exact invocation the human will type. This is the user-facing surface.

For each parameter, specify: **name, type, units, valid range, default, required?**

| Flag | Type | Required | Units | Range | Default | Description |
|------|------|----------|-------|-------|---------|-------------|
| -i --input | path | required | NA | NA | NA | Path to STEP file defining the lattice bounds |
| -o --output | path | optional | NA | NA | NA | Path and name of the output .step file, otherwise generated (e.g. input_file_name-cc5t1|
| -cc | float | required | mm | 0.4 - 50 | 5 | Distance between the bottom nodes of two adjecent cells |
| -t | float | required | mm | 0.4 - 20 | 1 | Side length of the diamond rod profile |
| -bg --background | flag | optional | NA | NA | disabled | Run worker processes at below-normal priority to reduce desktop impact |
| -v --verbose | flag | optional | NA | NA | disabled | Enable verbose console diagnostics while always writing a full `.log` file. |

**Exit:** 

Upon success the script shall produce an end of run summary report in the .log file and to console independent of the -verbose flag. This shall include: 
 - The runs input parameters 
 - Date and time of run start
 - Duration from start to completion
 - Run characteristics (number of tiles, number of parallel threads per stage of the generation procedure, etc)
 - Maximum memory usage
 - Path to output .step file 

Upon failure, the script shall output a human readable reason for the failure, e.g.: parameter bounds exceeded, issues with input geometry, issues with resarouces from the run-time system, write or read access issues, etc.
  
**Logging:**

A log file should be produced every run with the same name as the output file which is generated from the input file or provided by the -o flag. The log file should end with `.log` and should not include .step (that is only the last name for the geometry file.  

---

## 4. Geometry Domain Specification

### 4.1 Lattice unit cell type
The base geometry for the lattice is a strut-based uniform grid forming cube-like cells standing on its tip. The struts form the boundries of the cells along each edge. The struts have the profile of a square on its edge, like a diamond. The dimentions of the square profile is defined by input parameter `t`. Upon inspection of the end result, the struts are reclined from the normal axis (Z-axis) in degrees by the following calculation in numpy: np.degrees(np.arcsin(np.sqrt(2/3)))
Make sure to sure the native exact definition in native Jula language. It should be close to 55 degrees (but not exactly).
The distance between base points of each cell on the XY plane is defined by input parameter `cc`.
Upon inspection the rods protruding up from the xy plane from an intersecting node are separated by an angle of 120 degrees around the Z-axis.

### 4.2 Parametrization
- The sides of the diamond shaped square rod is defined as `t` in millimaters
- The diagonal distance between nodes across a cube is defined as `cc` in millimeters
- The bounds of the generated lattice and its placement in the xyz coordinate system is defined by the input step file. 

### 4.3 Boundary / shell requirements
- The lattice shall not have an outer solid shell generated around the build volume. It will be merged with the outer shell upon import into the enveloping part. However the truts must be closed against the geometry of provided input STEP file.
- No fillets/chamfers at strut junctions or bounding geometry

### 4.4 Performance & Optimization

Since this involves computational geometry:
- Profile geometry generation routines to identify bottlenecks
- Consider vectorization or parallelization for lattice tiling
- Cache expensive calculations (e.g., basis functions for triply-periodic surfaces)
- If caching to disk is used, put the files in a temporary folder `temp/<date><time>` where the output file is generated to. Clean up after a sucessful run. Leave for error analysis if the fun failes.

---

## 5. STEP Output Requirements
- **Clean up** Never produce a floating body smaller than`t` millimeters cubed.
  
- **STEP schema/AP:** AP203

- **Geometry representation in the file:** exact B-rep solid
  
- **Units:** mm
 
- **Metadata to embed:** Part name as concetenated <input_file_name>+cc<cc>+t<t> and generation parameters as STEP header

- **Downstream tool(s) that will open this file:** Soldiworks and Catia

---

## 6. Autonomous End-to-End Verification


### 6.1 Test scenarios
List concrete parameter sets that must be run automatically (at minimum: one small
case, one large/dense case, one edge case at parameter boundaries, one expected-failure
case for invalid input).

| Scenario | Parameters | Expected result |
|----------|-----------|------------------|
| smoke-fast | -i test/80mm-test-ball.step -cc 10 -t 2 | generation < 60 sec, quick test not applicable for output geometry verification |
| [TODO: for later] smoke-verified | -i test/80mm-test-ball.step -cc 10 -t 2 | valid STEP, generation < 60 sec, matching golden sample |
| [TODO: for later] dense-lattice | ... | valid STEP, no self-intersections, matching golden sample |
| [TODO: for later] invalid-input | wall-thickness > cell-size | exits nonzero, no file written |

### 6.2 Automated pass/fail checks
For every scenario the harness must verify, without human intervention:
- [ ] Process exits with expected console output
- [ ] STEP file is written and non-empty
- [ ] STEP file parses back successfully (round-trip read with the same or an
      independent library)
- [ ] Geometry is a valid closed manifold solid (no open edges / non-manifold edges)
- [ ] No self-intersections
- [ ] Bounding box of output matches requested `--input` within tolerance
- [ ] Runtime stays under an agreed performance budget for each scenario size (TBD until baseline is established)
- [ ] If golden sample is defined, check similarity of geometriers by running a geometry check script that inspects if the output geometry is inhibiting the same volume as the golden sample (e.g. subtraction either way should leave near zero volume). 


### 6.3 How verification runs offline
- Verification runs only in the dev/CI environment.
- Test runner: Test scripts and assets are separated into separate `tools/` folder.
- Results are reported as console summary for analysis and addition to the pull-request.

---

## 7. Error Handling & Edge Cases

- Invalid/out-of-range parameters should be rejected before any computation starts.
- If the 
- Read and write failures should be reported and result in a hard fail. Existing files can be overwritten.

---

## 8. Non-Functional Requirements

TBD

---

## 9. Open Questions / Decisions Needed

*Anything you're unsure about — list it here explicitly so it doesn't get silently
assumed by default. Delete each line once resolved.*

-
-


