# Conformal mapping of the human proximal femur for shape matching

<img src="https://github.com/NathanielNarra/Femur-RicciFlow/blob/master/docs/FemurIntroFig.png" align="right" width="450">

## Overview

This a repository for the MATLAB code developed for Ricci-flow based conformal paramterisation of the human proximal femur surface meshes for shape matching purposes. The procedure (illustrated [here](https://github.com/NathanielNarra/Femur-RicciFlow/blob/master/docs/Illustrated_procedural_flow.pdf)) is tuned specifically for this anatomy. It can be extended to other anatomical surfaces, however the current MATLAB script is only applicable to genus 0 surfaces. The scripts have been developed for the work presented in the manuscript (under revision): 

`Narra N, Abe S, Dimitrov V, Nikander R, Kouhia R, Sievänen H, Hyttinen J. “Ricci-flow based conformal mapping of the proximal femur to identify exercise loading effects”. (Under review)`


## Contents

|    **Folder**    | **Description of contents** |
|----------------|------------|
| Sample Femur Meshes | 20 meshes from the dataset used in the manuscript. These ‘.ply’ files contain triangular surface meshes (~50,000 faces). The surfaces represent only the proximal aspect of the femur with a boundary distal to the lesser trochanter process. |
| Subroutines      | Contains all functions accessed by the two main functions: STEP1_Parametrisation and STEP2_Template_Matching.  |
| Template | Contains the data structure containing the developed template mesh (3D and parametrised Annulus). Required for constructing iso-topological meshes of N surfaces through surface matching. |

| **Script** |  | **Description of function** |
|-------|--------|-----------|
| STEP1_Parametrisation.m | |iterates through N femur surfaces and parametrises to an annulus.   |
|| INPUT | list of file paths to mesh files of multiple (or single) surface mesh(es). Currently the script accepts ply format. read_ply.m (written by Gabriel Peyre) is provided. |
| |OUTPUT | The script does not explicitly output any data structure. However, all data structures are saved. It also saves a MATLAB ‘.fig’ file which renders the 3D surface mesh with the detected features (femoral head and greater trochanter) annotated. |
|STEP2_Template_Matching | | reads-in each data structure stored by STEP1 and matches a pre-defined template mesh (source mesh) to it (target mesh). Iterates through N surface meshes in the dataset. |
| |INPUT |folder id or folder path where the N parametrisation data structures are stored for N surface meshes. |
| |OUTPUT |returns a data structure with the modified template mesh and the positions of its nodes on the target meshes.  The position is represented as barycentric coordinates with respect to the triangular face it is incident on in the target mesh. This is termed as ‘natural representation’. This representation can be used to interpolate/attach the warped source mesh nodes with features from the target mesh. |

