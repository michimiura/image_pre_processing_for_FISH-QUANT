# image_pre_processing_for_FISH-QUANT
## 1. Run ImageJ macro Image_pre_processing_v2.0.ijm

This ImageJ macro takes raw z-stacked images acquired with Nikon Ti2-E and MetaMorph (multidimentional acquisition), and will

(1) stack up channels (e.g. DAPI, GFP, RFP and Cy5) from the same field of view,

(2) identify the focal plane and remove the extra z-stacks (focal plane +/- specified number of stacks),

(3) create a maximum projection from (2) and save,

(4) split the channels from (2) for FISH-QUANT and save each of the channels,

(5) create DAPI and GFP maximum projection of focal z-stacks (focal plane +/- 2) and save for segmentation.

## 2. Run Cellpose

For example,
```
python -m cellpose --dir /mnt/g/yyyymmdd/Segmentation --img_filt _DAPI --pretrained_model cyto --chan 0 --diameter 75 --save_tif --no_npy
python -m cellpose --dir /mnt/g/yyyymmdd/Segmentation --img_filt _GFP --pretrained_model cyto --chan 0 --diameter 150 --save_tif --no_npy
```
I find that the nuclei are more accurately captured with the option `--pretrained_model cyto` for my cells.

## 3. Run R script outline_from_cellpose.v3.R

For example from R Console,
```
source("outline_from_cellpose.v3.R")
```

This R script takes masks produced with cellpose (nuclei and cytoplasm), traces the perimeter of these masks with Moore's neighbour tracing algorithm, and produces a text file containing XY coordinates of the outlines (nucleus and cytoplasm) for each cell, in the format that is compatible with FISH-QUANT.  
