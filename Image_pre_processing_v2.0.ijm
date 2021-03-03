//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Image pre-processing written by Michi Miura, The University of Hong Kong, 09/03/2020
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// This macro takes raw z-stacked images acquired with Nikon Ti2-E and MetaMorph (multidimentional acquisition), and will
/// (1) stack up channels (e.g. DAPI, GFP, RFP and Cy5) from the same field of view,
/// (2) identify the focal plane and remove the extra z-stacks (focal plane +/- specified number of stacks),
/// (3) create a maximum projection from (2) and save,
/// (4) split the channels from (2) for FISH-QUANT and save each of the channels,
/// (5) create DAPI and GFP maximum projection of focal z-stacks (focal plane +/- 2) and save for segmentation
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//setBatchMode(true);

// select a folder
image_folder = getDirectory("Select an image folder to process");

// indicate the acquired channels
DAPI_w = getString("Indicate the channel for DAPI", "w1, w2, w3, w4");
GFP_w  = getString("Indicate the channel for GFP", "w1, w2, w3, w4");
RFP_w  = getString("Indicate the channel for RFP", "w1, w2, w3, w4 or NA");
Cy5_w  = getString("Indicate the channel for Cy5", "w1, w2, w3, w4 or NA");

// specify how many z-stacks +/- focal plane to keep.
z_lower = getString("Number of stacks to keep below focal plane", "18");
z_upper = getString("Number of stacks to keep above focal plane", "20");

// create folders
File.makeDirectory(image_folder + "Stack/");
File.makeDirectory(image_folder + "Stack/MAX/");
File.makeDirectory(image_folder + "FQ/");
File.makeDirectory(image_folder + "FQ/images/");

if(!(GFP_w=="NA")){File.makeDirectory(image_folder + "FQ/outlines_GFP/");}
if(!(RFP_w=="NA")){File.makeDirectory(image_folder + "FQ/outlines_RFP/");}
if(!(Cy5_w=="NA")){File.makeDirectory(image_folder + "FQ/outlines_Cy5/");}

if(!(GFP_w=="NA")){File.makeDirectory(image_folder + "FQ/results_GFP/");}
if(!(RFP_w=="NA")){File.makeDirectory(image_folder + "FQ/results_RFP/");}
if(!(Cy5_w=="NA")){File.makeDirectory(image_folder + "FQ/results_Cy5/");}

File.makeDirectory(image_folder + "Segmentation/");

// identify unique fields (DAPI)
f = getFileList(image_folder);
count = 0;
for (i=0; i<f.length; i++){
    
    if(endsWith(f[i], ".TIF")){ // select image files
        if(!matches(f[i], ".*thumb.*")) { // ignore thumb nails
            if(matches(f[i], ".*DAPI.*")) { // select DAPI images
                if(count==0){ // store the file name
                    field_id = f[i]; 
                    count++;
                }else{
                    field_id = Array.concat(field_id, f[i]); 
                    count++;
                }
            }
        }
    }
}

print("Processed with the ImageJ macro Image_pre_processing_v2.0.ijm");
print("Raw image name (DAPI)"+"\t"+"total slices"+"\t"+"focal plane"+"\t"+"lower slice"+"\t"+"upper slice"+"\t"+"New image name");


// look into each unique field
for (i=0; i < field_id.length; i++){

    flag = 0;

    // set the file name
    DAPI = field_id[i];
    if(!(GFP_w=="NA")){ GFP = replace( replace(field_id[i], "DAPI", "GFP") , DAPI_w, GFP_w ); }
    if(!(RFP_w=="NA")){ RFP = replace( replace(field_id[i], "DAPI", "RFP") , DAPI_w, RFP_w ); }
    if(!(Cy5_w=="NA")){ Cy5 = replace( replace(field_id[i], "DAPI", "Cy5") , DAPI_w, Cy5_w ); }

    // open the images
    open(image_folder + DAPI);
    if(!(GFP_w=="NA")){ open(image_folder + GFP); }
    if(!(RFP_w=="NA")){ open(image_folder + RFP); }
    if(!(Cy5_w=="NA")){ open(image_folder + Cy5); }
    
    // stack the channels
    if(!(GFP_w=="NA") && !(RFP_w=="NA") && !(Cy5_w=="NA")){       // GFP + RFP + Cy5
        run("Merge Channels...", "c1=[" + DAPI + "] c2=[" + GFP + "] c3=[" + RFP + "] c4=[" + Cy5 + "] create");
        Stack.setChannel(1); run("Cyan");
        Stack.setChannel(2); run("Green");
        Stack.setChannel(3); run("Yellow");
        Stack.setChannel(4); run("Magenta");
    }else if(!(GFP_w=="NA") && RFP_w=="NA" && !(Cy5_w=="NA")){    // GFP + Cy5
        run("Merge Channels...", "c1=[" + DAPI + "] c2=[" + GFP + "] c3=[" + Cy5 + "] create");
        Stack.setChannel(1); run("Cyan");
        Stack.setChannel(2); run("Green");
        Stack.setChannel(3); run("Magenta");
    }else if(!(GFP_w=="NA") && !(RFP_w=="NA") && Cy5_w=="NA"){    // GFP + RFP
        run("Merge Channels...", "c1=[" + DAPI + "] c2=[" + GFP + "] c3=[" + RFP + "] create");
        Stack.setChannel(1); run("Cyan");
        Stack.setChannel(2); run("Green");
        Stack.setChannel(3); run("Yellow");
    }
    //}else if(!(GFP_w=="NA") && RFP_w=="NA" && Cy5_w=="NA"){       // GFP
    //    run("Merge Channels...", "c1=[" + DAPI + "] c2=[" + GFP + "] create");
    //    Stack.setChannel(1); run("Cyan");
    //    Stack.setChannel(2); run("Green");
    //}else if(GFP_w=="NA" && !(RFP_w=="NA") && Cy5_w=="NA"){       // RFP
    //    run("Merge Channels...", "c1=[" + DAPI + "] c2=[" + RFP + "] create");
    //    Stack.setChannel(1); run("Cyan");
    //    Stack.setChannel(2); run("Yellow");
    //}else if(GFP_w=="NA" && RFP_w=="NA" && !(Cy5_w=="NA")){       // Cy5
    //    run("Merge Channels...", "c1=[" + DAPI + "] c2=[" + Cy5 + "] create");
    //    Stack.setChannel(1); run("Cyan");
    //    Stack.setChannel(2); run("Magenta");
    //}
    
    // the produced image should be "Composite"

    // extract dimensions
    Stack.getDimensions(width, height, channels, slices, frames); 

    for (z=1; z<=slices; z++){ // move through z-stacks and record the stdev of DAPI intensity in each frame

        Stack.setSlice(z);   // set the frame
        Stack.setChannel(1); // set to DAPI channel
        getStatistics(area, mean, min, max, std);

        if(z==1){x = std;}else{x = Array.concat(x, std);} // store the std in an object x

    }

    rankPosArr = Array.rankPositions(x); // an array of x coordinates, with x having been sorted from small to large 

    focal_plane = rankPosArr[rankPosArr.length -1] + 1; // slice with the largest std (slice index = x coordinate + 1)
    lower_boundary = focal_plane - z_lower; if(lower_boundary < 1){flag = 1;}
    upper_boundary = focal_plane + z_upper; if(upper_boundary > slices){flag = 1;}

    // adjust the contrast
    Stack.setSlice(focal_plane);
    for (chan=1; chan<=channels; chan++){
        Stack.setChannel(chan); run("Enhance Contrast", "saturated=0.35");
    }
    
    if(flag==0){

        // remove the excess z-stacks
        run("Duplicate...", "duplicate slices=lower_boundary-upper_boundary");
        
        // the produced image should be "Composite-1"

        // focal plane in the new stack
        new_focal_plane = focal_plane - lower_boundary + 1;

        // count the number of existing files with the same name
        current_list = getFileList(image_folder + "Stack/");
        test_name = replace( substring( DAPI, 0, lastIndexOf(DAPI, "_") ), "_" + DAPI_w + "WF DAPI Single", "" ); 
        test_name = substring(test_name, 0, lastIndexOf(test_name, "_"));

        suffix=1;
        for (hoi=0; hoi<current_list.length; hoi++){
            if(startsWith(current_list[hoi], test_name)){suffix++;}
        }

        // name the file with the suffix and save
        if(suffix < 10){name = test_name + "_0" + suffix;}
        if(suffix >= 10){name = test_name + "_" + suffix;}
        saveAs("tiff", image_folder + "Stack/" + name);

        // save the maximum projection
        run("Z Project...", "projection=[Max Intensity]");
        name_max_proj = getTitle();
        saveAs("Tiff", image_folder + "Stack/MAX/" + name_max_proj);

        // split the channels for FQ
        selectWindow(name + ".tif");
        run("Split Channels");
            
        // DAPI
        selectWindow("C1-" + name + ".tif"); 
        saveAs("tiff", image_folder + "FQ/images/" + name+"_DAPI");
        
        // for nuc segmentation
        run("Z Project...", "start=" + new_focal_plane - 2 + " stop=" + new_focal_plane + 2 +" projection=[Max Intensity]");
        name_DAPI_proj = replace(getTitle(), "MAX_", "ZPROJ_");
        saveAs("Tiff", image_folder + "Segmentation/" + name_DAPI_proj);
        
        // GFP
        selectWindow("C2-" + name + ".tif"); 
        saveAs("tiff", image_folder + "FQ/images/" + name+"_GFP");

        // for cyto segmentation
        run("Z Project...", "start=" + new_focal_plane -5 - 2 + " stop=" + new_focal_plane -5 + 2 +" projection=[Max Intensity]"); // offset 1.25 um (250 nm x 5)
        name_GFP_proj = replace(getTitle(), "MAX_", "ZPROJ_");
        saveAs("Tiff", image_folder + "Segmentation/" + name_GFP_proj);

        if(!(RFP_w=="NA") && !(Cy5_w=="NA")){
            // RFP
            selectWindow("C3-" + name + ".tif"); 
            saveAs("tiff", image_folder + "FQ/images/" + name+"_RFP");
            // Cy5
            selectWindow("C4-" + name + ".tif");
            saveAs("tiff", image_folder + "FQ/images/" + name+"_Cy5");
        }else if(RFP_w=="NA" && !(Cy5_w=="NA")){
            // Cy5
            selectWindow("C3-" + name + ".tif");
            saveAs("tiff", image_folder + "FQ/images/" + name+"_Cy5");
        }else if(!(RFP_w=="NA") && Cy5_w=="NA"){
            // RFP
            selectWindow("C3-" + name + ".tif"); 
            saveAs("tiff", image_folder + "FQ/images/" + name+"_RFP");
        }

        // close
        while(nImages>0){selectImage(nImages); close();} 
        run("Collect Garbage");

        print(image_folder + DAPI + "\t" + slices +"\t" + focal_plane + "\t" + lower_boundary + "\t" + upper_boundary + "\t" + image_folder + "Stack\\" + name +".tif");

    }

    if(flag==1){
        
        while(nImages>0){selectImage(nImages); close();} 
        run("Collect Garbage");
        print(image_folder + DAPI + "\t" + slices +"\t" + focal_plane + "\t" + lower_boundary + "\t" + upper_boundary + "\t" + "omitted");
        
    }

}

selectWindow("Log");
saveAs("txt", image_folder + "_Log");
close("Log");
