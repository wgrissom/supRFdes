# supRFdes
Code for implementation and validation of MRI RF pulse design by supervised learning

# Scripts
buildSARReg.m: Build SAR/VOP regularization matrix for a pulse design problem with a given set of VOP's and # of RF time points/slices
compare_kNN.m: The mother supervised learning shim comparison script. Runs all the methods, assuming msshim.m has already been run to generate the training shims.
loadData.m: Loads the maps and features. Invoked by compare_kNN.m.
msShim_randStart.m: Called by msshim.m to generate the training shims. Includes SAR regularization.
msShim_randStart_POCS_vopReg.m: Does shimming with POCS and VOP/SAR regularization. Called by msShim_POCS.m.
msshim.m: Loads b1+ maps+VOPs and runs msShim_randStart.m to get training shims. 
msshim_POCS.m: Calls msShim_randStart_POCS_vopReg.m to perform POCS shimming with VOP regularization. Is called by compare_kNN.m, but can also be run on its own for a given set of subjects/slices.
