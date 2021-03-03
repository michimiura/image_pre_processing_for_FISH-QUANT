my_dir <- readline(prompt="Enter folder name (yyyymmdd): ")

xy_calib <- readline(prompt="Enter scale (um/pixel): ")
xy_calib <- as.numeric(xy_calib)

primary_channel <- readline(prompt="Primary smFISH channel (RFP or Cy5): ")



### load libraries
library(EBImage)
library(ijtiff)
library(raster)
library(sp)
library(rgdal)

library(foreach)
library(doSNOW)
library(tcltk)

library(tidyverse)



### (def function) Moore neighbourhood for tracing boundary
moore <- function (m, len){

    #specify the number of cpu
    cl <- makeSOCKcluster(4)
    registerDoSNOW(cl)

    #set the progress bar
    # pb <- txtProgressBar(max=max(m), style=3)
    pb <- txtProgressBar(max=len, style=3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress=progress)
    
    # outlines <- vector("list", max(m))
    outlines <- vector("list", len)
    
    # outlines <- foreach(i = 1:max(m), .options.snow=opts)%dopar%{
    outlines <- foreach(i = 1:len, .options.snow=opts)%dopar%{

        region <- which(m==i, arr.ind=TRUE)
        if(length(region)==0){NULL}
        else{
        
        x <- min(region[,1])
        y <- max(region[which(region[,1]==x),2])
    
        X <- x
        Y <- y
    
        trans <- "N"
    
        flag <- 0
        while(flag == 0){
    
            if(trans == "N"){
                if     (m[x-1, y+1] == i){ X <- c(X, x-1); Y <- c(Y, y+1); x <- x-1; y <- y+1; trans <- "SW" }
                else if(m[x-1, y  ] == i){ X <- c(X, x-1); Y <- c(Y, y);   x <- x-1; y <- y;   trans <- "W"  }
                else if(m[x-1, y-1] == i){ X <- c(X, x-1); Y <- c(Y, y-1); x <- x-1; y <- y-1; trans <- "NW" }
                else if(m[x,   y-1] == i){ X <- c(X, x);   Y <- c(Y, y-1); x <- x;   y <- y-1; trans <- "N"  }
                else if(m[x+1, y-1] == i){ X <- c(X, x+1); Y <- c(Y, y-1); x <- x+1; y <- y-1; trans <- "NE" }
                else if(m[x+1, y  ] == i){ X <- c(X, x+1); Y <- c(Y, y);   x <- x+1; y <- y;   trans <- "E"  }
                else if(m[x+1, y+1] == i){ X <- c(X, x+1); Y <- c(Y, y+1); x <- x+1; y <- y+1; trans <- "SE" }
                else                     { X <- c(X, x);   Y <- c(Y, y+1); x <- x;   y <- y+1; trans <- "S"  }
            }else if(trans == "NE"){
                if     (m[x-1, y  ] == i){ X <- c(X, x-1); Y <- c(Y, y);   x <- x-1; y <- y;   trans <- "W"  }
                else if(m[x-1, y-1] == i){ X <- c(X, x-1); Y <- c(Y, y-1); x <- x-1; y <- y-1; trans <- "NW" }
                else if(m[x,   y-1] == i){ X <- c(X, x);   Y <- c(Y, y-1); x <- x;   y <- y-1; trans <- "N"  }
                else if(m[x+1, y-1] == i){ X <- c(X, x+1); Y <- c(Y, y-1); x <- x+1; y <- y-1; trans <- "NE" }
                else if(m[x+1, y  ] == i){ X <- c(X, x+1); Y <- c(Y, y);   x <- x+1; y <- y;   trans <- "E"  }
                else if(m[x+1, y+1] == i){ X <- c(X, x+1); Y <- c(Y, y+1); x <- x+1; y <- y+1; trans <- "SE" }
                else if(m[x,   y+1] == i){ X <- c(X, x);   Y <- c(Y, y+1); x <- x;   y <- y+1; trans <- "S"  }
                else                     { X <- c(X, x-1); Y <- c(Y, y+1); x <- x-1; y <- y+1; trans <- "SW" }
            }else if(trans == "E"){
                if     (m[x-1, y-1] == i){ X <- c(X, x-1); Y <- c(Y, y-1); x <- x-1; y <- y-1; trans <- "NW" }
                else if(m[x,   y-1] == i){ X <- c(X, x);   Y <- c(Y, y-1); x <- x;   y <- y-1; trans <- "N"  }
                else if(m[x+1, y-1] == i){ X <- c(X, x+1); Y <- c(Y, y-1); x <- x+1; y <- y-1; trans <- "NE" }
                else if(m[x+1, y  ] == i){ X <- c(X, x+1); Y <- c(Y, y);   x <- x+1; y <- y;   trans <- "E"  }
                else if(m[x+1, y+1] == i){ X <- c(X, x+1); Y <- c(Y, y+1); x <- x+1; y <- y+1; trans <- "SE" }
                else if(m[x,   y+1] == i){ X <- c(X, x);   Y <- c(Y, y+1); x <- x;   y <- y+1; trans <- "S"  }
                else if(m[x-1, y+1] == i){ X <- c(X, x-1); Y <- c(Y, y+1); x <- x-1; y <- y+1; trans <- "SW" }
                else                     { X <- c(X, x-1); Y <- c(Y, y);   x <- x-1; y <- y;   trans <- "W"  }
            }else if(trans == "SE"){
                if     (m[x,   y-1] == i){ X <- c(X, x);   Y <- c(Y, y-1); x <- x;   y <- y-1; trans <- "N"  }
                else if(m[x+1, y-1] == i){ X <- c(X, x+1); Y <- c(Y, y-1); x <- x+1; y <- y-1; trans <- "NE" }
                else if(m[x+1, y  ] == i){ X <- c(X, x+1); Y <- c(Y, y);   x <- x+1; y <- y;   trans <- "E"  }
                else if(m[x+1, y+1] == i){ X <- c(X, x+1); Y <- c(Y, y+1); x <- x+1; y <- y+1; trans <- "SE" }
                else if(m[x,   y+1] == i){ X <- c(X, x);   Y <- c(Y, y+1); x <- x;   y <- y+1; trans <- "S"  }
                else if(m[x-1, y+1] == i){ X <- c(X, x-1); Y <- c(Y, y+1); x <- x-1; y <- y+1; trans <- "SW" }
                else if(m[x-1, y  ] == i){ X <- c(X, x-1); Y <- c(Y, y);   x <- x-1; y <- y;   trans <- "W"  }
                else                     { X <- c(X, x-1); Y <- c(Y, y-1); x <- x-1; y <- y-1; trans <- "NW" }
            }else if(trans == "S"){
                if     (m[x+1, y-1] == i){ X <- c(X, x+1); Y <- c(Y, y-1); x <- x+1; y <- y-1; trans <- "NE" }
                else if(m[x+1, y  ] == i){ X <- c(X, x+1); Y <- c(Y, y);   x <- x+1; y <- y;   trans <- "E"  }
                else if(m[x+1, y+1] == i){ X <- c(X, x+1); Y <- c(Y, y+1); x <- x+1; y <- y+1; trans <- "SE" }
                else if(m[x,   y+1] == i){ X <- c(X, x);   Y <- c(Y, y+1); x <- x;   y <- y+1; trans <- "S"  }
                else if(m[x-1, y+1] == i){ X <- c(X, x-1); Y <- c(Y, y+1); x <- x-1; y <- y+1; trans <- "SW" }
                else if(m[x-1, y  ] == i){ X <- c(X, x-1); Y <- c(Y, y);   x <- x-1; y <- y;   trans <- "W"  }
                else if(m[x-1, y-1] == i){ X <- c(X, x-1); Y <- c(Y, y-1); x <- x-1; y <- y-1; trans <- "NW" }
                else                     { X <- c(X, x);   Y <- c(Y, y-1); x <- x;   y <- y-1; trans <- "N"  }
            }else if(trans == "SW"){
                if     (m[x+1, y  ] == i){ X <- c(X, x+1); Y <- c(Y, y);   x <- x+1; y <- y;   trans <- "E"  }
                else if(m[x+1, y+1] == i){ X <- c(X, x+1); Y <- c(Y, y+1); x <- x+1; y <- y+1; trans <- "SE" }
                else if(m[x,   y+1] == i){ X <- c(X, x);   Y <- c(Y, y+1); x <- x; y <- y+1;   trans <- "S"  }
                else if(m[x-1, y+1] == i){ X <- c(X, x-1); Y <- c(Y, y+1); x <- x-1; y <- y+1; trans <- "SW" }
                else if(m[x-1, y  ] == i){ X <- c(X, x-1); Y <- c(Y, y);   x <- x-1; y <- y;   trans <- "W"  }
                else if(m[x-1, y-1] == i){ X <- c(X, x-1); Y <- c(Y, y-1); x <- x-1; y <- y-1; trans <- "NW" }
                else if(m[x,   y-1] == i){ X <- c(X, x);   Y <- c(Y, y-1); x <- x; y <- y-1;   trans <- "N"  }
                else                     { X <- c(X, x+1); Y <- c(Y, y-1); x <- x+1; y <- y-1; trans <- "NE" }
            }else if(trans == "W"){
                if     (m[x+1, y+1] == i){ X <- c(X, x+1); Y <- c(Y, y+1); x <- x+1; y <- y+1; trans <- "SE" }
                else if(m[x,   y+1] == i){ X <- c(X, x);   Y <- c(Y, y+1); x <- x;   y <- y+1; trans <- "S"  }
                else if(m[x-1, y+1] == i){ X <- c(X, x-1); Y <- c(Y, y+1); x <- x-1; y <- y+1; trans <- "SW" }
                else if(m[x-1, y  ] == i){ X <- c(X, x-1); Y <- c(Y, y);   x <- x-1; y <- y;   trans <- "W"  }
                else if(m[x-1, y-1] == i){ X <- c(X, x-1); Y <- c(Y, y-1); x <- x-1; y <- y-1; trans <- "NW" }
                else if(m[x,   y-1] == i){ X <- c(X, x);   Y <- c(Y, y-1); x <- x;   y <- y-1; trans <- "N"  }
                else if(m[x+1, y-1] == i){ X <- c(X, x+1); Y <- c(Y, y-1); x <- x+1; y <- y-1; trans <- "NE" }
                else                     { X <- c(X, x+1); Y <- c(Y, y);   x <- x+1; y <- y;   trans <- "E"  }
            }else if(trans == "NW"){
                if     (m[x,   y+1] == i){ X <- c(X, x);   Y <- c(Y, y+1); x <- x;   y <- y+1; trans <- "S"  }
                else if(m[x-1, y+1] == i){ X <- c(X, x-1); Y <- c(Y, y+1); x <- x-1; y <- y+1; trans <- "SW" }
                else if(m[x-1, y  ] == i){ X <- c(X, x-1); Y <- c(Y, y);   x <- x-1; y <- y;   trans <- "W"  }
                else if(m[x-1, y-1] == i){ X <- c(X, x-1); Y <- c(Y, y-1); x <- x-1; y <- y-1; trans <- "NW" }
                else if(m[x,   y-1] == i){ X <- c(X, x);   Y <- c(Y, y-1); x <- x;   y <- y-1; trans <- "N"  }
                else if(m[x+1, y-1] == i){ X <- c(X, x+1); Y <- c(Y, y-1); x <- x+1; y <- y-1; trans <- "NE" }
                else if(m[x+1, y  ] == i){ X <- c(X, x+1); Y <- c(Y, y);   x <- x+1; y <- y;   trans <- "E"  }
                else                     { X <- c(X, x+1); Y <- c(Y, y+1); x <- x+1; y <- y+1; trans <- "SE" }
            }
        
            if(X[1] == X[length(X)] && Y[1] == Y[length(Y)]){flag <- 1}

        }
        
        list(X, Y)
        
        }
    
    }
    close(pb)
    stopCluster(cl)
    
    outlines

}
### (end of function)



### (def function) add a scale bar
add_scale <- function(l){
    plot(c((1508-l/xy_calib),1508), c(100,100), xlim=c(1,L), ylim=c(1,L), type="l", lend = "butt", col="white", lwd = 2,
         xaxs = "i", yaxs = "i", xlab=NA, ylab=NA, xaxt="n", yaxt="n")
}
### (end of function)



### (def function) find the centre of mass for each mask
find_centre <- function(in_mask){
    
    cell_ind <- as.numeric(in_mask) %>% unique() %>% sort() %>% tail(n = -1)
    
    if(!length(cell_ind) >= 1){NULL}
    else{
        
        for(i in 1:length(cell_ind)){
        
            mass <- sum(in_mask == cell_ind[i])
            
            colSum <- colSums(in_mask == cell_ind[i])
            # colSum_cum <- 0
            # for (iii in 1:dim(mask_nuc)[1]){
            # colSum_cum <- colSum_cum + colSum[iii]
            # if(colSum_cum >= mass/2){nuc_c.x <- iii; break}
            # }
            if(i == 1){nuc.x <- sum(colSum * 1:L) / mass}
            else{nuc.x <- c(nuc.x, sum(colSum * 1:L) / mass)}

            rowSum <- rowSums(in_mask == cell_ind[i])
            # rowSum_cum <- 0
            # for (iii in 1:dim(mask_nuc)[2]){
            # rowSum_cum <- rowSum_cum + rowSum[iii]
            # if(rowSum_cum >= mass/2){nuc_c.y <- iii; break}
            # }
            if(i == 1){nuc.y <- sum(rowSum * 1:L) / mass}
            else{nuc.y <- c(nuc.y, sum(rowSum * 1:L) / mass)}
            
        }
        
        data.frame(idx = cell_ind, x = nuc.x, y = nuc.y)
        
    }
    
}
### (end of function)



###
f_DAPI <- list.files(path = paste0("G:/", my_dir, "/Segmentation/"), pattern = "_DAPI_cp_masks.tif", full.names = TRUE)

for(i in 1:length(f_DAPI)){
#for(i in 23){
    
    flag <- 0
    
    ### read the mask (nuc)
    mask_nuc  <- read_tif(f_DAPI[i])[,,1,1] %>% as.Image() %>% EBImage::transpose()
    storage.mode(mask_nuc) <- "integer"
    
    ### read the mask (cyto)
    mask_cyto <- read_tif(gsub("DAPI", "GFP", f_DAPI[i]))[,,1,1] %>% as.Image() %>% EBImage::transpose()
    storage.mode(mask_cyto) <- "integer"
    
    ### read images
    img_DAPI <- read_tif(gsub("_cp_masks", "", f_DAPI[i])) %>% as_EBImage()
    img_cyto <- read_tif(gsub("_cp_masks", "", gsub("DAPI", "GFP", f_DAPI[i])))[,,1,1] %>% as_EBImage()
    
    ### width
    L <- dim(mask_nuc)[1]
    
    ### display the microscopy image
    mic_img <- rgbImage(red = img_DAPI * 1/quantile(img_DAPI, 0.8),
                        green = img_cyto * 1/quantile(img_cyto, 0.8), 
                        blue  = img_DAPI * 1/quantile(img_DAPI, 0.8))
    #graphics.off(); windows(width=17/2.54, height=17/2.54); par(mfrow=c(4,4)); par(mar=c(0,0,0,0))
    graphics.off(); windows(width=10/2.54, height=20/2.54); par(mfrow=c(4,2)); par(mar=c(0,0,0,0))
    display(mic_img)
    par(new=TRUE)
    add_scale(20) # scale bar = 20 um
    
    ###
    plot(0,0, xaxs = "i", yaxs = "i", xlab=NA, ylab=NA, xaxt="n", yaxt="n", col = "white")
    #plot(0,0, xaxs = "i", yaxs = "i", xlab=NA, ylab=NA, xaxt="n", yaxt="n")
    #plot(0,0, xaxs = "i", yaxs = "i", xlab=NA, ylab=NA, xaxt="n", yaxt="n")
    
    ### display original cyto masks
    display(colorLabels(mask_cyto))
    
    ### add cyto index
    c.mass <- find_centre(mask_cyto)
    par(new=TRUE)
    plot(c.mass$y, L - c.mass$x, xlim=c(1,L), ylim=c(1, L), #note the tricky coordinate conversion
         xaxs = "i", yaxs = "i", xlab=NA, ylab=NA, xaxt="n", yaxt="n",
         pch = 16, col = "red")
    text(c.mass$y, L - c.mass$x, c.mass$idx, cex = 1, pos = 3,col = "red")
    
    ### remove cytoplasm on the edge
    ind_rm <- c(mask_cyto[1,],     #left
                mask_cyto[L,],     #right
                mask_cyto[,1],     #top
                mask_cyto[,L]) %>% #bottom
              unique() %>%
              sort()
    mask_cyto <- rmObjects(mask_cyto, ind_rm)
    
    ### display updated cyto mask
    display(colorLabels(mask_cyto))
    
    ### add cyto index
    c.mass <- find_centre(mask_cyto)
    par(new=TRUE)
    plot(c.mass$y, L - c.mass$x, xlim=c(1,L), ylim=c(1, L), #note the tricky coordinate conversion
         xaxs = "i", yaxs = "i", xlab=NA, ylab=NA, xaxt="n", yaxt="n",
         pch = 16, col = "red")
    text(c.mass$y, L - c.mass$x, c.mass$idx, cex = 1, pos = 3,col = "red")
    
    ### display original nuc mask
    display(colorLabels(mask_nuc))
    
    ### add nuc index
    c.mass <- find_centre(mask_nuc)
    par(new=TRUE)
    plot(c.mass$y, L - c.mass$x, xlim=c(1,L), ylim=c(1, L), #note the tricky coordinate conversion
         xaxs = "i", yaxs = "i", xlab=NA, ylab=NA, xaxt="n", yaxt="n",
         pch = 16, col = "red")
    text(c.mass$y, L - c.mass$x, c.mass$idx, cex = 1, pos = 3,col = "red")
    
    ### remove nucleus on the edge
    ind_rm <- c(mask_nuc[1,],     #left
                mask_nuc[L,],     #right
                mask_nuc[,1],     #top
                mask_nuc[,L]) %>% #bottom
              unique() %>%
              sort()
    mask_nuc <- rmObjects(mask_nuc, ind_rm)
    
    ### display updated nuc mask
    display(colorLabels(mask_nuc))
    
    ### add nuc index
    c.mass <- find_centre(mask_nuc)
    par(new=TRUE)
    plot(c.mass$y, L - c.mass$x, xlim=c(1,L), ylim=c(1, L), #note the tricky coordinate conversion
         xaxs = "i", yaxs = "i", xlab=NA, ylab=NA, xaxt="n", yaxt="n",
         pch = 16, col = "red")
    text(c.mass$y, L - c.mass$x, c.mass$idx, cex = 1, pos = 3,col = "red")
    
    ###
    gc()
    
    if(!is.na(c.mass)){
    
        ### assign cyto index to nuc masks
        c.mass <- cbind(c.mass, cyto_idx = mask_cyto[cbind(c.mass$y, c.mass$x)])
        for(ii in 1:nrow(c.mass)){
        
            mask_nuc[which(mask_nuc == c.mass[ii,]$idx)] <- c.mass[ii,]$cyto_idx * 1000
            
        }
        mask_nuc <- mask_nuc/1000
        
        ### new centre of mass
        c.mass.new <- find_centre(mask_nuc)
        
        if(!is.na(c.mass.new)){
        
            ### display
            paintObjects(mask_cyto, 
                         paintObjects(mask_nuc, mic_img, col=c(rgb(1,1,1),NA), thick=TRUE, closed=TRUE),
                         col=c(rgb(1,1,1),NA), thick=TRUE, closed=TRUE) %>% 
            display()
            par(new=TRUE); par(mar=c(0,0,0,0))
            plot(c.mass.new$y, L - c.mass.new$x, xlim=c(1,L), ylim=c(1, L), #note the tricky coordinate conversion
                 xaxs = "i", yaxs = "i", xlab=NA, ylab=NA, xaxt="n", yaxt="n",
                 pch = 16, col = "red")
            text(c.mass.new$y, L - c.mass.new$x, c.mass.new$idx, cex = 1, pos = 3,col = "red")
            
            ###
            gc()
        
            ### correct nuclei if it goes beyond the cytoplasm
            for(ii in 1:nrow(c.mass.new)){
        
                tmp_nuc <- which(mask_nuc == c.mass.new[ii,]$idx)
                tmp_cyto_out <- which(!mask_cyto == c.mass.new[ii,]$idx)
                ind_rm <- intersect(tmp_nuc, tmp_cyto_out)
                if(!length(ind_rm) == 0){mask_nuc[ind_rm] <- 0}
            
            }
            
            ### produce an outline text file for FQ (Moore)
            my_max <- max(max(mask_nuc), max(mask_cyto))
            print("Tracing the nuc boundary..."); flush.console()
            outlines_nuc <- moore(mask_nuc, my_max)
            print("Tracing the cyto boundary..."); flush.console()
            outlines_cyto <- moore(mask_cyto, my_max)
            
            flag <- 1
    
            ### check the traces
            graphics.off(); windows(width=8.5/2.54, height=8.5/2.54); par(mar=c(0,0,0,0)); par(pty="s")
            display(mic_img)
            par(new=TRUE)
            add_scale(20) # scale bar = 20 um
            for (ii in 1:length(outlines_nuc)){

                par(new=TRUE)
                plot(outlines_nuc[[ii]][[1]], outlines_nuc[[ii]][[2]], xlim=c(1,L), ylim=c(L,1),
                     type="l", xaxs="i", yaxs="i", xaxt="n", yaxt="n", xlab=NA, ylab=NA, col="blue")
                par(new=TRUE)
                plot(outlines_cyto[[ii]][[1]], outlines_cyto[[ii]][[2]], xlim=c(1,L), ylim=c(L,1),
                     type="l", xaxs="i", yaxs="i", xaxt="n", yaxt="n", xlab=NA, ylab=NA, col="green")
                
            }

            ### save
            rp <- vector("list", 1); dev.set(2); rp[[1]] <- recordPlot()
            tiff(paste0("G:/", my_dir, "/Segmentation/", gsub("_DAPI_cp_masks.tif", "_outline.tif", basename(f_DAPI[i]))),
                 width = 8.5/2.54, height = 8.5/2.54, units = "in", res = 300)
            replayPlot(rp[[1]])
            dev.off()
            
        }
        
    }
    
    ### write outline.txt
    if(primary_channel == "Cy5"){
        outline_name <- paste0("G:/", my_dir, "/FQ/outlines_Cy5/",
                               gsub(  "_DAPI_cp_masks.tif", "", gsub("ZPROJ_", "", basename(f_DAPI[i]))  ),
                               "_Cy5_outline.txt")
    }else if(primary_channel == "RFP"){
        outline_name <- paste0("G:/", my_dir, "/FQ/outlines_RFP/",
                               gsub(  "_DAPI_cp_masks.tif", "", gsub("ZPROJ_", "", basename(f_DAPI[i]))  ),
                               "_RFP_outline.txt")
    }
    
    sink(outline_name)

    ### imitate the outline descriptions created by FQ and CellProfiler (e.g.)
    
    # FISH-QUANT	v3a
    # File-version	3D_v1
    # RESULTS OF SPOT DETECTION PERFORMED ON 17-Feb-2020 
    # COMMENT	
    # IMG_Raw	C4-200212_S4_Composite_s1.tif
    # IMG_Filtered	
    # IMG_DAPI	C1-200212_S4_Composite_s1.tif
    # IMG_TS_label	
    # FILE_settings	
    # PARAMETERS
    # Pix-XY	Pix-Z	RI	Ex	Em	NA	Type
    # 110	300	1.515	647	670	1.45	widefield
    # CELL_START	Cell_2
    # X_POS	1136	1138	1143	1149	1152	1159	1160	1163	1165	1170	1173	1177	1180	1182	1187	1192	1195	1198	1201	1203	1208	1211	1214	1216	1217	1219	1222	1225	1227	1228	1233	1237	1239	1240	1240	1241	1241	1244	1246	1246	1247	1248	1251	1252	1252	1253	1254	1255	1256	1256	1256	1258	1258	1258	1258	1258	1258	1258	1258	1258	1258	1258	1258	1258	1258	1258	1258	1258	1258	1258	1258	1258	1258	1258	1258	1258	1258	1258	1259	1258	1258	1257	1256	1256	1255	1253	1250	1247	1244	1242	1241	1241	1239	1235	1227	1222	1219	1216	1214	1214	1210	1205	1201	1199	1197	1192	1189	1180	1174	1167	1160	1149	1146	1143	1141	1135	1131	1125	1121	1115	1110	1107	1105	1102	1096	1093	1090	1088	1085	1079	1078	1076	1074	1072	1071	1068	1063	1060	1058	1054	1050	1048	1042	1036	1033	1032	1031	1029	1027	1027	1022	1016	1012	1008	1005	1001	997	996	996	995	991	988	985	981	979	977	975	971	970	967	966	963	962	962	960	957	953	951	951	951	949	949	949	949	949	949	949	949	949	949	949	949	950	951	951	952	952	953	955	955	955	956	957	958	958	958	959	960	963	964	966	966	969	970	972	974	977	981	983	986	990	993	996	997	1002	1006	1010	1012	1016	1018	1022	1025	1028	1033	1036	1041	1043	1047	1053	1057	1063	1067	1070	1073	1083	1091	1100	1103	1104	1107	1111	1113	1115	1117	1118	1120	1124	1126	1127	1128	1129	1131	1132	1134	1134	
    # Y_POS	782	782	783	784	784	784	784	785	785	787	788	789	792	796	799	800	801	803	805	808	811	814	815	817	817	817	820	822	825	826	831	834	837	840	842	842	844	847	851	855	857	860	864	866	868	870	873	874	876	878	882	885	886	888	891	898	902	904	906	909	912	917	920	922	924	929	932	937	942	947	951	955	966	967	970	976	979	983	991	995	996	998	1000	1001	1003	1007	1011	1014	1018	1020	1022	1023	1026	1029	1035	1039	1042	1045	1047	1048	1054	1055	1058	1059	1059	1062	1066	1069	1072	1074	1076	1079	1080	1080	1082	1084	1085	1087	1089	1091	1093	1095	1096	1098	1099	1100	1100	1101	1102	1104	1104	1105	1105	1105	1105	1105	1105	1107	1107	1107	1108	1108	1108	1108	1108	1108	1108	1108	1108	1109	1109	1109	1109	1109	1109	1107	1107	1106	1106	1105	1102	1099	1096	1092	1090	1088	1085	1081	1077	1075	1074	1071	1071	1070	1067	1059	1050	1046	1044	1040	1034	1031	1029	1028	1024	1019	1012	1005	998	992	986	981	974	969	966	961	954	951	945	937	932	928	922	917	912	910	907	899	894	888	882	879	875	872	868	862	857	852	849	845	838	833	829	825	821	817	814	812	809	808	806	805	804	801	800	800	800	799	797	795	794	793	792	792	790	787	785	784	784	784	783	783	783	783	783	783	783	783	783	783	783	783	783	783	783	
    # Z_POS	
    # CELL_END

    cat("FISH-QUANT\tv3a\n")
    cat("File-version\t3D_v1\n")
    cat("RESULTS OF SPOT DETECTION PERFORMED ON 15-Aug-1945\n")
    cat("COMMENT\t\n")
    if(primary_channel == "Cy5"){
        cat(paste("IMG_Raw\t", gsub("_DAPI_cp_masks", "_Cy5", gsub("ZPROJ_", "", basename(f_DAPI[i]))), "\n", sep=""))
    }else if(primary_channel == "RFP"){
        cat(paste("IMG_Raw\t", gsub("_DAPI_cp_masks", "_RFP", gsub("ZPROJ_", "", basename(f_DAPI[i]))), "\n", sep=""))
    }
    cat("IMG_Filtered\t\n")
    cat(paste("IMG_DAPI\t", gsub("_cp_masks", "", gsub("ZPROJ_", "", basename(f_DAPI[i]))), "\n", sep=""))
    cat(paste("IMG_TS_label\t", gsub(  "DAPI", "GFP", gsub("_cp_masks", "", gsub("ZPROJ_", "", basename(f_DAPI[i])))  ), "\n", sep=""))
    cat("FILE_settings\t\n")
    cat("PARAMETERS\n")
    cat("Pix-XY\tPix-Z\tRI\tEx\tEm\tNA\tType\n")
    cat("110\t300\t1.515\t647\t670\t1.45\twidefield\n")
    
    if(flag == 1){
        
        # for(ii in 1:length(outlines_nuc)){
        for(ii in 1:my_max){

            cat(  paste("CELL_START\t", "Cell_", ii, "\n", sep="")  )
            cat(  paste("X_POS\t", paste(  head(outlines_cyto[[ii]][[1]],-1), collapse="\t"  ), "\t\n", sep="")  )
            cat(  paste("Y_POS\t", paste(  head(outlines_cyto[[ii]][[2]],-1), collapse="\t"  ), "\t\n", sep="")  )
            cat("Z_POS\t\n")
            cat("CELL_END\n")
            cat("Nucleus_START\tNuc_CP\n")
            cat(  paste("X_POS\t", paste(  head(outlines_nuc[[ii]][[1]],-1), collapse="\t"  ), "\t\n", sep="")  )
            cat(  paste("Y_POS\t", paste(  head(outlines_nuc[[ii]][[2]],-1), collapse="\t"  ), "\t\n", sep="")  )
            cat("Z_POS\t\n")
            cat("Nucleus_END\n")

        }
        
    }

    sink()
    
    if(i==length(f_DAPI)){print("Done!"); flush.console()}
    
    ###
    gc()

}