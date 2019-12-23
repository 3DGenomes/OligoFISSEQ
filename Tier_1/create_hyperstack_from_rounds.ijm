//===============================================================================
// Creates a single hyperstack from rounds folder and correct 3D drift.
//===============================================================================
macro "Create hyperstack from rounds"  {
 dirPath = getDirectory( "Choose a Directory" ); 
 FolderTemplatePath=dirPath+"round_";
 TotoPath=dirPath+"toto";
 outputPath = dirPath+"merged"+File.separator;
 seg_channel=5;
 toto_seg_channel=5;
 y_drift_c2=-1;
 create_hyperstack_from_rounds(FolderTemplatePath, TotoPath, outputPath, seg_channel, toto_seg_channel, y_drift_c2);
}

function create_hyperstack_from_rounds(FolderTemplatePath, TotoPath, outputPath, seg_channel, toto_seg_channel, y_drift_c2) {
  //setBatchMode(true);
  round_nbr=1;
  list_round1 = getFileList(FolderTemplatePath+round_nbr);
  File.makeDirectory(outputPath);
  for (i = 0; i < list_round1.length; i++) {
    numImage=nImages();
  //for (i = 0; i < 2; i++) {
    max_z=0;
    round_nbr=1;
    while(File.exists(FolderTemplatePath+round_nbr+File.separator+list_round1[i])) {
    //for (p = 0; p < FolderListPath.length; p++) {
      namest=toLowerCase(list_round1[i]);
      if (endsWith(namest,".tif") || endsWith(namest,".tiff")) { 
		open(FolderTemplatePath+round_nbr+File.separator+list_round1[i]);
      } else {
		run("Bio-Formats", "open='"+FolderTemplatePath+round_nbr+File.separator+list_round1[i]+"'"+
		"autoscale color_mode=Grayscale view=Hyperstack");
      }
      print('Check maximum z');
      print(FolderTemplatePath+round_nbr+File.separator+list_round1[i]);
      getDimensions(width2, height2, channelCount2, sliceCount2, frameCount2);
      if (sliceCount2>max_z) max_z = sliceCount2;
      close();
      round_nbr+=1;
    }

	//toto
    filenameToto = File.nameWithoutExtension;
    find_toto = false;
    if(File.exists(TotoPath+File.separator+filenameToto+".tiff")) {
      open(TotoPath+File.separator+filenameToto+".tiff");
      print(TotoPath+File.separator+filenameToto+".tiff");
      find_toto = true;
    } else if(File.exists(TotoPath+File.separator+filenameToto+".tif")) {
      open(TotoPath+File.separator+filenameToto+".tif");
      print(TotoPath+File.separator+filenameToto+".tif");
      find_toto = true;
    } else if(File.exists(TotoPath+File.separator+filenameToto+".nd2")) {
      run("Bio-Formats", "open='"+TotoPath+File.separator+filenameToto+".nd2"+"'"+
      "autoscale color_mode=Grayscale view=Hyperstack");
      print(TotoPath+File.separator+filenameToto+".nd2");
      find_toto = true;
    }

    if(find_toto) {
      getDimensions(width2, height2, channelCount2, sliceCount2, frameCount2);
      if (sliceCount2>max_z) max_z = sliceCount2;
      close();
    }
    
    zslices=max_z;
    print('Maximum z:'+zslices);
    stack_str = "";
    round_nbr=1;
    while(File.exists(FolderTemplatePath+round_nbr+File.separator+list_round1[i])) {
      namest=toLowerCase(list_round1[i]);
      if (endsWith(namest,".tif") || endsWith(namest,".tiff")) { 
		open(FolderTemplatePath+round_nbr+File.separator+list_round1[i]);
      } else {
		run("Bio-Formats", "open='"+FolderTemplatePath+round_nbr+File.separator+list_round1[i]+"'"+
		"autoscale color_mode=Grayscale view=Hyperstack");
      }
      print(FolderTemplatePath+round_nbr+File.separator+list_round1[i]);
      selectImage(numImage+round_nbr);
      title2 = getTitle();
      filename = getFilename();
      getDimensions(width2, height2, channelCount2, sliceCount2, frameCount2);
      
      // correct channel2 chromatic shift
      print('Correcting channel 2 chromatic shift');
      for (d = 0; d < sliceCount2*channelCount2; d+=channelCount2) {
		setSlice(d+2);
		run("Translate...", "x=0 y="+y_drift_c2+" interpolation=None slice");
      }
      
      //match dimensions in z
      for (d = 0; d < (zslices-sliceCount2); d++) {
		duptitle = "M" + sliceCount2 + "_" + title2;
		run("Duplicate...", "duplicate slices="+sliceCount2+" title='"+duptitle+"'");
		run("Concatenate...", "  image1='" + title2 + "' image2='" + duptitle + "' title='"+title2+"'");
      }
      selectImage(numImage+round_nbr);
      rename(filename+"_round"+round_nbr);
      stack_str = stack_str + " stack"+round_nbr+"='"+filename+"_round"+round_nbr+"'";
	  round_nbr+=1;
    }

    //toto
    filenameToto = File.nameWithoutExtension;
    find_toto = false;
    if(File.exists(TotoPath+File.separator+filenameToto+".tiff")) {
      open(TotoPath+File.separator+filenameToto+".tiff");
      print(TotoPath+File.separator+filenameToto+".tiff");
      find_toto = true;
    } else if(File.exists(TotoPath+File.separator+filenameToto+".tif")) {
      open(TotoPath+File.separator+filenameToto+".tif");
      print(TotoPath+File.separator+filenameToto+".tif");
      find_toto = true;
    } else if(File.exists(TotoPath+File.separator+filenameToto+".nd2")) {
      run("Bio-Formats", "open='"+TotoPath+File.separator+filenameToto+".nd2"+"'"+
      "autoscale color_mode=Grayscale view=Hyperstack");
      print(TotoPath+File.separator+filenameToto+".nd2");
      find_toto = true;
    }
    
    if(find_toto) {
      print('Adding toto image');
      getDimensions(width2, height2, channelCount_toto, sliceCount_toto, frameCount_toto);
      //match dimensions in z
      for (d = 0; d < (zslices-sliceCount_toto); d++) {
		duptitle = "M" + sliceCount_toto + "_" + title2;
		run("Duplicate...", "duplicate slices="+sliceCount_toto+" title='"+duptitle+"'");
		run("Concatenate...", "  image1='" + title2 + "' image2='" + duptitle + "' title='"+title2+"'");
      }
      //match dimensions in c
      concstr = "";
      for (d = 0; d < channelCount2; d++) {
      	if(d < (channelCount2-channelCount_toto))
			concstr = concstr + "c"+(d+1)+"=[C1-" + title2 + "] ";
		else 
			concstr = concstr + "c"+(d+1)+"=[C"+(d-(channelCount2-channelCount_toto)+1)+"-"+ title2 + "] ";
      }
      print(concstr);
	  run("Split Channels");
	  run("Merge Channels...", concstr+" create");
      
      selectImage(numImage+round_nbr);
      rename(filename+"_toto");
      stack_str = stack_str + " stack"+round_nbr+"='"+filename+"_toto'";
	  round_nbr+=1;
    }
    round_nbr-=1;
    selectImage(filename+"_round1");
    getVoxelSize(pxwidth, pxheight, pxdepth, unit);
    run("Concatenate...", stack_str);
    
    // there is only one image left; rename it to the original title
    rename(filename);
    selectImage(filename);
    
    // set the correct dimensions
    Stack.setDimensions(channelCount2, zslices, round_nbr);
    selectImage(filename);
    print("Detecting nuclei to compute 3D drift");
    run("Z Project...", "projection=[Max Intensity]");
    setOption("BlackBackground", true);
    run("Convert to Mask", "method=Default background=Dark calculate black");
    setSlice(seg_channel);
    run("Median...", "radius=2 slice");
    run("Find Edges", "slice");
    run("Analyze Particles...", "size=10.00-Infinity add slice");
    close();
    selectImage(filename);
    //selectImage(numImage+2);
    roiManager("Combine");
    //run("Correct 3D drift", "channel=4 multi_time_scale sub_pixel edge_enhance only=0 lowest=1 highest="+zslices);
    print('Using channel '+seg_channel+' to correct 3D drift');
    run("Correct 3D drift", "channel="+seg_channel+" only=0 lowest=1 highest="+zslices);
    
    run("Select None");
    roiManager("Delete");
    selectImage(filename);
    close();
    //selectImage(numImage+2);
    selectImage("registered time points");
    rename(filename);
    getDimensions(width2, height2, channelCount2, sliceCount2, frameCount2);
    zslices=sliceCount2;
    //remove empty slices
    slicesperround=channelCount2*zslices;
    max_start_offset = 0;
    for (r = 0; r < round_nbr; r++) {
      start_offset=1;
      for (d = 0; d < zslices; d++) {
	setSlice(d+channelCount2*d+slicesperround*r+1);
	getStatistics(n, mean, min, max, std, histogram);
	if(max==0) start_offset = d + 1;
	else break;
      }
      if(start_offset>max_start_offset) max_start_offset=start_offset;
    }
    max_end_offset = 0;
    for (r = 0; r < round_nbr; r++) {
      end_offset=0;
      for (d = 0; d < zslices; d++) {
	//print(channelCount2*(zslices-d)+slicesperround*r);
	setSlice(channelCount2*(zslices-d)+slicesperround*r);
	getStatistics(n, mean, min, max, std, histogram);
	if(max==0) end_offset=d+1;
	else break;
      }
      if(end_offset>max_end_offset) max_end_offset=end_offset;
      //print(max_end_offset);
    }
    max_end_offset = zslices - max_end_offset;
    print('Start non-empty image:'+max_start_offset);
    print('End non_empty image:'+max_end_offset);
    selectImage(filename);
    //run("Make Substack...", "channels=1-"+channelCount2+" slices="+(max_start_offset+1)+"-"+max_end_offset+" frames=1-"+FolderListPath.length);
    run("Make Substack...", "channels=1-"+channelCount2+" slices="+(max_start_offset+1)+"-"+max_end_offset+" frames=1-"+round_nbr);
    run("Properties...", "unit="+unit+" pixel_width="+pxwidth+" pixel_height="+pxheight+" voxel_depth="+pxdepth);

	//saveAs("tiff", outputPath + filename);
    rename(filename+"_full");
    full_image = getImageID();
    selectImage(filename);
    close();
    selectImage(full_image);
    getDimensions(width2, height2, channelCount2, sliceCount2, frameCount2);
	File.makeDirectory(outputPath + filename);

    // save roi
    //run("Duplicate...", "duplicate");
    //proj_image = getImageID();	  
    run("Z Project...", "projection=[Max Intensity]");
    setOption("BlackBackground", true);
    run("Convert to Mask", "method=Default background=Dark calculate black");
    setSlice(seg_channel);
    run("Median...", "radius=2 slice");
    run("Watershed", "slice");
    //run("Find Edges", "slice");
    run("Analyze Particles...", "size=100.00-Infinity add slice");
    //roiManager("Save", outputPath + filename + "_roi" + ".zip");
    close();
    run("Select None");
    
	stack_str_full="";
	for (r = 0; r < frameCount2; r++) {
		stack_str="";
		for (c = 0; c < channelCount2; c++) {
		  selectImage(full_image);
		  run("Make Substack...", "channels="+(c+1)+" slices="+1+"-"+sliceCount2+" frames="+(r+1));
		  rename(filename+"_full_"+(r+1)+"_"+(c+1));
	      run("Hyperstack to Stack");
	      selectImage(filename+"_full_"+(r+1)+"_"+(c+1));
		  run("Enhance Contrast...", "saturated=0 equalize process_all use");
		  stack_str = stack_str + " image"+(c+1)+"='"+filename+"_full_"+(r+1)+"_"+(c+1)+"'";
	    }
	    run("Concatenate...", stack_str);
	    run("Stack to Hyperstack...", "order=xyzct channels="+channelCount2+" slices="+sliceCount2+" frames=1 display=Grayscale");
	    //Stack.setDimensions(channelCount2, sliceCount2, 1);
	    rename(filename+"_full_"+(r+1));
	    stack_str_full = stack_str_full + " image"+(r+1)+"='"+filename+"_full_"+(r+1)+"'";
	}
	//print(stack_str_full);
	run("Concatenate...", stack_str_full);
	Stack.setDimensions(channelCount2, sliceCount2, frameCount2);

	//getVoxelSize(pixelWidth, pixelHeight, pixelDepth, cal_unit);
	
	for (u=0; u<roiManager("count"); ++u) {
		//Roi.getBounds(rx, ry, rwidth, rheight);
		run("Duplicate...", "duplicate title=crop");
        roiManager("Select", u);
        run("To Bounding Box");
        run("Enlarge...", "enlarge=8 pixel");
        run("Crop");
	    roiManager("Select", u);
	    Roi.move(8, 8);
	    run("Enlarge...", "enlarge=1"); //sometimes we cut too much, we increase 1 
	    roiManager("Update");
		roiManager("save selected", outputPath + filename + File.separator + 'roi_'+(u+1)+'_roi.zip');
	    run("Select None");
	    run("Remove Overlay");
	    setVoxelSize(pxwidth, pxheight, pxdepth, unit);
        run("Properties...", "unit="+unit+" pixel_width="+pxwidth+" pixel_height="+pxheight+" voxel_depth="+pxdepth);
        saveAs("tiff", outputPath + filename + File.separator + 'roi_'+(u+1));
        close();
        selectImage(filename+"_full");
    }
    //roiManager("Save", outputPath + filename + "_roi" + ".zip");
    roiManager("Delete");
    //close();
    //selectImage(filename);
    //close();
     while (nImages>numImage) { 
      selectImage(nImages); 
      close(); 
     }
  }
  print('Done!'); 
}


macro "=========================="{} 

function getFilename() {
  t = getTitle();
  ext = newArray(".tif", ".tiff", ".lif", ".lsm", ".czi", ".nd2", ".ND2");    
  for(i=0; i<ext.length; i++)
    t = replace(t, ext[i], "");  
  return t;
}