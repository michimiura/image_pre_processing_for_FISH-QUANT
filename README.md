# image_pre_processing
/// This macro takes raw z-stacked images acquired with Nikon Ti2-E and MetaMorph (multidimentional acquisition), and will
/// (1) stack up channels (e.g. DAPI, GFP, RFP and Cy5) from the same field of view,
/// (2) identify the focal plane and remove the extra z-stacks (focal plane +/- specified number of stacks),
/// (3) create a maximum projection from (2) and save,
/// (4) split the channels from (2) for FISH-QUANT and save each of channel,
/// (5) create DAPI and GFP maximum projection of selected z-stacks (focal plane +/- 2) and save for segmentation.
