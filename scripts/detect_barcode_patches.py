from __future__ import division
from os import path, listdir, walk

from ij import IJ, ImagePlus, ImageStack, CompositeImage

from fiji.plugin.trackmate import Model
from fiji.plugin.trackmate import Settings
from fiji.plugin.trackmate import TrackMate
from fiji.plugin.trackmate import Logger

from fiji.plugin.trackmate.detection import LogDetectorFactory
from fiji.plugin.trackmate.detection import DogDetectorFactory
from fiji.plugin.trackmate.tracking.sparselap import SparseLAPTrackerFactory
from fiji.plugin.trackmate.tracking import LAPUtils
from fiji.plugin.trackmate.util import TMUtils

import fiji.plugin.trackmate.Dimension as Dimension
import fiji.plugin.trackmate.features.FeatureFilter as FeatureFilter
import fiji.plugin.trackmate.features.spot.SpotIntensityAnalyzerFactory as SpotIntensityAnalyzerFactory
import fiji.plugin.trackmate.features.track.TrackDurationAnalyzer as TrackDurationAnalyzer
import fiji.plugin.trackmate.io.TmXmlWriter as TmXmlWriter

import ij.plugin.Duplicator as Duplicator
import ij.plugin.ContrastEnhancer as ContrastEnhancer
from ij.process import AutoThresholder
from ij.process.AutoThresholder import getThreshold
from ij.process.AutoThresholder import Method
from ij.plugin.filter import ParticleAnalyzer as PA
from ij.plugin.filter import EDM as EDM
from ij.plugin.frame import RoiManager
from ij.measure import ResultsTable
from ij.plugin import ZProjector
from ij.process import ShortProcessor as ShortProcessor
from ij.process import ByteProcessor as ByteProcessor

import org.jfree.chart.renderer.InterpolatePaintScale as InterpolatePaintScale
import java.awt.Color as Color
import java.io.File as File
import sys
import re
import csv
import math
import random
from jarray import zeros
from itertools import tee

#@ File(label="Select a directory containing input images", style="directory") input_folder
#@ File(label="Select a barcode file") barcode_file
#@ Integer(label="Number of rounds",value=5) nbr_rounds
#@ String(label="Subsample on rounds",value=1234) rounds_sub
#@ Integer(label="Number of channels",value=4) nchannels
#@ Integer(label="Channel used to segment nuclei if no roi file provided",value=5) segment_channel
#@ Boolean(label="Use toto",value=0) use_toto
#@ Integer(label="Channel used in toto image",value=4) toto_channel

max_displacement_xy = 0.3 # 0.267 microns/pixel?
max_displacement_z = 0.6
max_displacement_xy_merge = 0.6 # 0.267 microns/pixel?
min_nuclei_size = 350 # pixel^2
max_nuclei_size = 100000 # pixel^2
nbr_rounds_sub = len(rounds_sub)

radius_puncta = 0.5
show_barcoding=False
show_intensity_image=True
convert_channel = {0:'0',1:'4',2:'3',3:'2',4:'1',5:'0'}
max_float_value = (2**16)/2-1
barcode_models = []

colorToRGB = {
'Red' : [255,0,0],
'Green' : [0,255,0],
'Blue' : [0,0,255],
'Orange' : [255,127,0],
'Cyan' : [0,255,255],
'Yellow' : [255,255,0],
'Magenta' : [255,0,255],
'Indigo' : [75,0,130],
'Violet' : [238,130,238],
'Greyscale' : [255,255,255],
'Aquamarine' : [127,255,212],
'Navy Blue' : [0,0,128],
'Sky Blye' : [135,206,235],
'Turquoise' : [64,224,208],
'Beige' : [245,245,220],
'Brown' : [165,42,42],
'Chocolate' : [210,105,30],
'Dark wood' : [133,94,66],
'Light wood' : [133,99,99],
'Olive' : [128,128,0],
'Green yellow' : [173,255,47],
'Sea green' : [32,178,170],
'Khaki' : [240,230,140],
'Salmon' : [250,128,114],
'Pink' : [255,192,203],
'Tomato' : [255,99,71],
'Scarlet' : [140,23,23],
'Purple' : [128,0,128],
'Wheat' : [245,222,179],
'Silver grey' : [192,192,192],
'Cerise': [222, 49, 99],
'Candy':[255,8,0],
'Fuchsia':[255,0,255],
'Bubble gum': [255, 193, 204],
'Berry': [255, 67, 164],
'Darkred': [139,0,0],
'FireBrick': [178, 34, 34],
'Crimson':[220,20,60],
'IndianRed':[205,92,92],
'DarkGreen': [0,153,0],
'LightBlue': [51,255,255],
'LightPurple': [204,204,255],
'LightBrown': [204,204,0],
'White': [255,255,255]
}

# Encode color RGB in floats:
tmp = {}
for c,rgb in colorToRGB.iteritems():
  tmp[c] = [v/255.0 for v in rgb]
colorToRGB = tmp

# Colors in the desired listing order:
colors = ['Red', 'Green', 'Blue',
		  'Purple', 'Brown', 'Yellow',
          'Orange', 'Silver grey',
          'Cyan', 'Magenta', 'Indigo',
          'Turquoise', 'Tomato', 'Olive',
          'Violet', 'Green yellow', 'Khaki',
          'Scarlet', 'Beige', 'Chocolate',
          'Pink', 'Wheat',
          'Sea green', 'Greyscale', 'Light wood',
          'Sky Blye',  'Salmon', 'Navy Blue',
          'Aquamarine', 'Dark wood', 'White']

chr_colors = [('DarkGreen','Green'),('Orange','Yellow'),('Blue','LightBlue'),
			  ('Purple','LightPurple'),('Darkred','Red'),('Greyscale','Silver grey')]
chr_mono_colors = ('Red','White')
# colors are floating-point 3-element tuples (doing RGBA is a matter of setting the component count)
def poly_gradient(colors, steps, components=3):

    def linear_gradient(start, finish, substeps):
        yield start
        for i in range(1, substeps):
            yield tuple([(start[j]+(float(i)/(substeps-1))*(finish[j]-start[j])) for j in range(components)])

    def pairs(seq):
        a, b = tee(seq)
        next(b, None)
        return zip(a, b)

    substeps = int(float(steps)/(len(colors)-1))

    for a, b in pairs(colors):
        for c in linear_gradient(a, b, substeps):
            yield c

def histogram(data, nBins):
	max_data, min_data = max(data), min(data)
	binWidth = (max_data - min_data) / nBins
	hist=[0 for x in xrange(nBins)]
	for i in xrange(len(data)):
		index=min( int(math.floor( ( data[ i ] - min_data ) / binWidth )), nBins - 1 )
		hist[index] += 1
	return hist

def center_point(p1, p2):

    return (int((p1[0] + p2[0])/2),int((p1[1] + p2[1])/2),int((p1[2] + p2[2])/2))

def adjacent(pixs,p1):

	return any([(((abs(p1[0]-p2[0])<=1 and p1[1]-p2[1] == 0) or (abs(p1[1]-p2[1])<=1 and p1[0]-p2[0] == 0)) and abs(p1[2]-p2[2])==0) \
				or (((abs(p1[1]-p2[1])==0 and p1[0]-p2[0] == 0)) and abs(p1[2]-p2[2])<=1) for p2 in pixs])

def median(lst):
    sortedLst = sorted(lst)
    lstLen = len(lst)
    index = (lstLen - 1) // 2

    if (lstLen % 2):
        return sortedLst[index]
    else:
        return (sortedLst[index] + sortedLst[index + 1])/2.0

def getDistance(a, b):

    if len(a) != len(b):
        raise Exception("ERROR: non comparable points")

    accumulatedDifference = 0.0
    for i in range(len(a)):
        squareDifference = pow((a[i]-b[i]), 2)
        accumulatedDifference += squareDifference

    return math.sqrt(accumulatedDifference)

def reduce_patch(patch_bar, zslices, max_area):
	max_offset = int(math.sqrt(max_area)/2.)
	pixs = patch_barcode['pixels']
	center_pix = patch_barcode['max_intensity_pixel']
	reduced_pixs = [pix for i in xrange(zslices) for pix in pixs \
	        if (abs(pix[0]-center_pix[0]) <= max_offset and (pix[1]-center_pix[1]) <= max_offset) and (pix[2] == i)]
	return reduced_pixs

def load_probes(input_file):
    probes = []
    with open(input_file, 'r') as csvfile:
        csvreader = csv.reader(csvfile, delimiter='\t')
        header = next(csvreader, None)
        for row in csvreader:
            probe = {}
            for i,h in enumerate(header):
                probe[h] = row[i]
            probes.append(probe)
    return probes

probes=load_probes(str(barcode_file))
barcodes = [int(probe['Readout']) for probe in probes]
locis = [probe['Probe'] for probe in probes]
chrs = []
loci_chroms = {}
for loci in locis:
    if 'p' in loci or 'q' in loci:
        chrom = (loci.split('p')[0]).split('q')[0]
    else:
        chrom = loci[0]
    if chrom in loci_chroms:
        loci_chroms[chrom].append(loci)
    else:
        loci_chroms[chrom] = [loci]
        chrs.append(chrom)

max_spots_in_nucleus = len(locis)*4

input_folder = str(input_folder)
for root, dirs, files in walk(input_folder):
    for file in files:
    	if file.endswith(".tiff") or file.endswith(".tif"):
			tiff_curr = path.join(root, file)
			print tiff_curr
			file_pattern, file_extension_pattern = path.splitext(tiff_curr)
			label_tiff=tiff_curr
			label_file_pattern=file_pattern
			imp = IJ.openImage(tiff_curr)
			nslices = imp.getNSlices()
			nframes = nbr_rounds
			width = imp.width
			height = imp.height
			pix_red = zeros(width * height, 'h')
			pix_red[0] = int(max_float_value)
			pixelWidth = imp.getCalibration().pixelWidth
			pixelDepth = imp.getCalibration().pixelDepth
			# Get normalization parameters
			p99 = {}
			pmode = {}
			max_global_intensity = 0
			stack = imp.getStack()
			for frame in range(1,nframes+1):
				for c in range(1,nchannels+1):
					pixels = []
					max_mode = 0
					max_max = 0
					for i in range(1,nslices+1):
						index = imp.getStackIndex(c, i, frame)
						stackProc = stack.getProcessor(index)
						stats = stackProc.getStatistics()
						if stats.max > max_max:
							max_max = stats.max
						if stats.mode > max_mode:
							max_mode = stats.mode
						pixels += stackProc.convertToFloat().getPixels()
					p99[(frame,c)] = max_max
					pmode[(frame,c)] = 0
					max_glob = 1.0
					if max_glob > max_global_intensity:
						max_global_intensity = max_glob

			roimanager = RoiManager(True)
			roiFile = '%s_roi.zip' % (file_pattern)
			if path.isfile(roiFile):
				print 'Using existing roi file %s' % roiFile
				roimanager.runCommand("Open",roiFile)
			else:
				# Limit the detection to nuclei
				channel1sliceImp = Duplicator().run(imp, segment_channel, segment_channel, 1, nslices, 1, 1)
				project = ZProjector()
				project.setMethod(ZProjector.MAX_METHOD)
				project.setImage(channel1sliceImp)
				project.doProjection()
				projection = project.getProjection()
				projection.setTitle('Projection '+imp.getTitle())
				IJ.run(projection, "Auto Threshold", "method=Huang ignore_black white")
				IJ.run(projection, "Convert to Mask", "")
				improc = projection.getProcessor()
				improc.medianFilter()
				dm = EDM()
				dm.toWatershed(improc)
				edges = ImagePlus('Edges', improc)
				edges.show()

				rt = ResultsTable()
				p = PA(PA.ADD_TO_MANAGER + PA.SHOW_OUTLINES, PA.AREA + PA.CENTER_OF_MASS, rt, min_nuclei_size, max_nuclei_size)
				p.setRoiManager(roimanager)
				p.setHideOutputImage(True)
				p.analyze(edges)
				mmap = p.getOutputImage()
				IJ.save(mmap, input_folder+'Nuclei '+imp.getTitle())
				edges.close()
			roi_count = roimanager.getCount()
			pixel_data = []
			for roi_idx in xrange(roi_count):
				pixel_data.append({})
				roi = roimanager.getRoi(roi_idx)
				points_roi = roi.getContainedPoints()
				for frame in xrange(nbr_rounds):
					print 'Compiling round %i\n' % (frame)
					# get max_intensity channel

					for i in xrange(nslices):
						print 'Compiling slice %i\n' % (i)
						for c in range(1,nchannels+1):
							index = imp.getStackIndex(c, i+1, frame+1)
							stackProc = stack.getProcessor(index).convertToFloat()
							pix = dict(((point_roi.x,point_roi.y),stackProc.getPixelValue(point_roi.x,point_roi.y)) for point_roi in points_roi)
							#pix = dict(((w,h),stackProc.getPixelValue(w,h)) for w in xrange(width) for h in xrange(height))
							for p in pix:
								if (p[0],p[1]) not in pixel_data[roi_idx]:
									#pixel_data[roi_idx][p] = [x[:] for x in [[0]*(nframes*2+1)]*nslices]
									pixel_data[roi_idx][p] = [{k:v for (k,v) in zip(['%s_%d'%(mi,fr)
						  for mi,fr in zip(['max_intensity',
								    'max_intensity_channel',
								    'second_intensity']*2*nframes,
								  [d for d in xrange(nframes) for _ in xrange(3)])],
						[0]*3*nframes)} for _ in xrange(nslices)]
								pix_intensity = 0  if p99[(frame+1,c)] == 0 \
									else (pix[p] - pmode[(frame+1,c)])/(p99[(frame+1,c)] - pmode[(frame+1,c)])

								if pixel_data[roi_idx][p][i]['second_intensity_%d'%frame] < pix_intensity:
									if pixel_data[roi_idx][p][i]['max_intensity_%d'%frame] <= pix_intensity:
										pixel_data[roi_idx][p][i]['second_intensity_%d'%frame] = pixel_data[roi_idx][p][i]['max_intensity_%d'%frame]
										pixel_data[roi_idx][p][i]['max_intensity_%d'%frame] = pix_intensity
										pixel_data[roi_idx][p][i]['max_intensity_channel_%d'%frame] = c
									else:
										pixel_data[roi_idx][p][i]['second_intensity_%d'%frame] = pix_intensity

			if use_toto:
				# toto
				totoImp = Duplicator().run(imp, toto_channel, toto_channel, 1, nslices, nbr_rounds+1, nbr_rounds+1)
				totoImp.copyScale(imp)
				IJ.save(totoImp, '%s_toto.tif' % (file_pattern))
				totoImpstack = totoImp.getStack()
				with open('%s_toto_intensities.tsv' % (file_pattern), mode='w') as toto_file:
				    tsv_writer = csv.writer(toto_file, delimiter='\t', quotechar='"', quoting=csv.QUOTE_MINIMAL)
				    tsv_writer.writerow(['roi','pixel_x','pixel_y','pixel_z','intensity'])
				    for roi_idx in xrange(roi_count):
						pixel_data.append({})
						roi = roimanager.getRoi(roi_idx)
						points_roi = roi.getContainedPoints()
						for i in xrange(nslices):
							totoImpProc = totoImpstack.getProcessor(i+1).convertToFloat()
							pix = dict(((point_roi.x,point_roi.y),totoImpProc.getPixelValue(point_roi.x,point_roi.y)) for point_roi in points_roi)
							for p in pix:
								tsv_writer.writerow([roi_idx,p[0],p[1],i+1,pix[p]])

			# create intensity image
			newstack = ImageStack(width, height)
			for i in xrange(nslices):
				pix_int = zeros(width * height, 'h')
				tmp_pix_int = [0 for _ in xrange(width * height)]
				for c in range(1,nchannels+1):
					index = imp.getStackIndex(c, i+1, 1)
					stackProc = stack.getProcessor(index)
					pixels = stackProc.convertToFloat().getPixels()
					pixels_norm = map(lambda x: 0  if p99[(1,c)] == 0 else (x - pmode[(1,c)])/(p99[(1,c)] - pmode[(1,c)]), pixels)
					tmp_pix_int = [tmp_pix_int[pix_i] if tmp_pix_int[pix_i] > pixels_norm[pix_i] else pixels_norm[pix_i] for pix_i in xrange(width * height)]
				pix_int = map(lambda x: int((x/max_global_intensity)*max_float_value), tmp_pix_int)
				newstack.addSlice(None, ShortProcessor(width, height, pix_int, None))

			imp_int = ImagePlus(imp.title+"_max_intensity", newstack)
			imp_int.setDimensions(1, nslices, 1)
			imp_int.copyScale(imp)
			IJ.save(imp_int, '%s_max_intensity.tif' % (file_pattern))

			# create barcode image
			stack_color = ImageStack(width, height)
			stack_int = ImageStack(width, height)
			pix_int = {}
			patches = {}
			patches[nbr_rounds] = [dict((barc,[]) for barc in barcodes) for roi_idx in xrange(roi_count)]
			patches[nbr_rounds_sub] = [dict((barc,[]) for barc in barcodes) for roi_idx in xrange(roi_count)]
			max_area_pix = 5
			base_colors_multi = [(colorToRGB[col_i],colorToRGB[col_e]) for col_i,col_e in chr_colors]
			base_colors_mono = (colorToRGB[chr_mono_colors[0]],colorToRGB[chr_mono_colors[1]])
			barcode_colors = {}
			for idc, chrom in enumerate(chrs):
			    base_colors = base_colors_mono if len(loci_chroms[chrom]) > 10 else base_colors_multi[idc]
			    grad_colors = [pol for pol in poly_gradient([base_colors[0],base_colors[1]],
			                                                       len(loci_chroms[chrom]))]
			    for loci_id, loci in enumerate(loci_chroms[chrom]):
			        barcode_colors[loci] = grad_colors[loci_id]

			for i in xrange(nslices):
				print 'Barcodind slices %i\n' % (i)
				pix_red = zeros(width * height, 'h')
				pix_green = zeros(width * height, 'h')
				pix_blue = zeros(width * height, 'h')
				for id_b, barcode in enumerate(barcodes):
					pix_int[barcode] = zeros(width * height, 'h')
				for roi_idx in xrange(roi_count):
				#for roi_idx in [0]:
					print 'Barcoding roi %i\n' % (roi_idx)
					for p in pixel_data[roi_idx]:
						#new_barcode = ''.join([convert_channel[pixel_data[roi_idx][p][i][frame]] for frame in xrange(nframes)])
						new_barcode = ''.join([convert_channel[pixel_data[roi_idx][p][i]['max_intensity_channel_%d'%frame]] for frame in xrange(nframes)])
						#pix_intensity = sum([(pixel_data[roi_idx][p][i][nframes+frame]/max_global_intensity) for frame in xrange(nframes)])/nframes
						pix_intensity = sum([(pixel_data[roi_idx][p][i]['max_intensity_%d'%frame]/max_global_intensity) for frame in xrange(nframes)])/nframes
						#pixel_data[roi_idx][p][i][2*nframes] = new_barcode
						pixel_data[roi_idx][p][i]['barcode'] = new_barcode
						for id_b, barcode in enumerate(barcodes):
						#for id_b, barcode in enumerate(tmp_barcode):
							match_barcodes = 0
							subs_barcode = ''.join([str(barcode)[int(barc_i)-1] for barc_i in str(rounds_sub)])
							subs_new_barcode = ''.join([str(new_barcode)[int(barc_i)-1] for barc_i in str(rounds_sub)])
							if subs_barcode == subs_new_barcode:
								match_barcodes = nbr_rounds_sub
							if str(barcode)[:nbr_rounds] == new_barcode[:nbr_rounds]:
								match_barcodes = nbr_rounds
							#match_barcodes = sum([a==b for a, b in zip(str(barcode)[:nbr_rounds],new_barcode[:nbr_rounds])])
							if match_barcodes >= nbr_rounds_sub :
								#print id_b, barcode, int(id_b/6), int(id_b%6)

								found_patches = []
								for id_p, patch_barcode in enumerate(patches[match_barcodes][roi_idx][barcode]):
									if adjacent(patch_barcode['pixels'],(p[0],p[1],i)):
										found_patches.append(id_p)

								if len(found_patches)>0:
									patch_barcode = patches[match_barcodes][roi_idx][barcode][found_patches[0]]
									#add pixel to patch
									patch_barcode['pixels'].append((p[0],p[1],i))
									if patch_barcode['max_intensity'] < pix_intensity:
										patch_barcode['max_intensity_pixel'] = (p[0],p[1],i)
										patch_barcode['max_intensity'] = pix_intensity
										for frame in xrange(nframes):
											patch_barcode['max_intensity_%d'%frame] = pixel_data[roi_idx][p][i]['max_intensity_%d'%frame]
											patch_barcode['second_intensity_%d'%frame] = pixel_data[roi_idx][p][i]['second_intensity_%d'%frame]

									#merge other patches
									for id_p in found_patches[1:]:
										patch_barcode['pixels'] +=  patches[match_barcodes][roi_idx][barcode][id_p]['pixels']
										if patch_barcode['max_intensity'] < patches[match_barcodes][roi_idx][barcode][id_p]['max_intensity']:
											patch_barcode['max_intensity_pixel'] = patches[match_barcodes][roi_idx][barcode][id_p]['max_intensity_pixel']
											patch_barcode['max_intensity'] = patches[match_barcodes][roi_idx][barcode][id_p]['max_intensity']
											for frame in xrange(nframes):
												patch_barcode['max_intensity_%d'%frame] = patches[match_barcodes][roi_idx][barcode][id_p]['max_intensity_%d'%frame]
												patch_barcode['second_intensity_%d'%frame] = patches[match_barcodes][roi_idx][barcode][id_p]['second_intensity_%d'%frame]

									patches[match_barcodes][roi_idx][barcode] = [pat for n_id_p,pat in enumerate(patches[match_barcodes][roi_idx][barcode]) if n_id_p not in found_patches[1:]]
								else:

									new_pix_dict = {'pixels': [(p[0],p[1],i)],
													'max_intensity': pix_intensity,
													'max_intensity_pixel': (p[0],p[1],i),
													'roi':roi_idx,
													'match_barcodes':match_barcodes}
									for frame in xrange(nframes):
										new_pix_dict['max_intensity_%d'%frame] = pixel_data[roi_idx][p][i]['max_intensity_%d'%frame]
										new_pix_dict['second_intensity_%d'%frame] = pixel_data[roi_idx][p][i]['second_intensity_%d'%frame]
									patches[match_barcodes][roi_idx][barcode].append(new_pix_dict)

								if match_barcodes == nbr_rounds:
									pix_int[barcode][p[1] * width + p[0]] = int(pix_intensity*max_float_value)
									grad_col = barcode_colors[locis[id_b]]
									pix_red[p[1] * width + p[0]] = int((max_float_value)*grad_col[0])
									pix_green[p[1] * width + p[0]] = int((max_float_value)*grad_col[1])
									pix_blue[p[1] * width + p[0]] = int((max_float_value)*grad_col[2])
				for id_b, barcode in enumerate(barcodes):
					stack_int.addSlice(None, ShortProcessor(width, height, pix_int[barcode], None))
				stack_color.addSlice(None, ShortProcessor(width, height, pix_red, None))
				stack_color.addSlice(None, ShortProcessor(width, height, pix_green, None))
				stack_color.addSlice(None, ShortProcessor(width, height, pix_blue, None))

			tmp_imp = ImagePlus(imp.title+"_barcoded", stack_color)
			tmp_imp.setDimensions(3, nslices, 1)
			comp_imp = CompositeImage(tmp_imp, CompositeImage.COMPOSITE)
			comp_imp.setCalibration(imp.getCalibration())
			IJ.save(comp_imp, '%s_barcoded.tif' % (file_pattern))
			if show_barcoding:
				comp_imp.show()
			#print pix_int[32311]
			tmp_imp = ImagePlus(imp.title+"_max_intensity", stack_int)
			tmp_imp.setDimensions(len(barcodes), nslices, 1)

			imp_int = CompositeImage(tmp_imp, CompositeImage.GRAYSCALE)
			imp_int.setCalibration(imp.getCalibration())
			IJ.save(imp_int, '%s_barcodes_splitted.tif' % (file_pattern))
			if show_barcoding:
				imp_int.show()

			#clean noise. big patches get reduced to max 40 pixels per plane around max intensity
			for match_barcodes in [nbr_rounds,nbr_rounds_sub]:
				for roi_idx in xrange(roi_count):
					for barcode in barcodes:
						for patch_barcode in patches[match_barcodes][roi_idx][barcode]:
							patch_barcode['pixels'] = reduce_patch(patch_barcode, nslices, 50)

			#export to csv
			with open('%s_punctas.tsv' % (file_pattern), mode='w') as punctas_file:
			    tsv_writer = csv.writer(punctas_file, delimiter='\t', quotechar='"', quoting=csv.QUOTE_MINIMAL)
			    tsv_writer.writerow(['barcode','roi','max_intensity_pixel_x','max_intensity_pixel_y',
								  'max_intensity_pixel_z','max_intensity', 'nbr_pixels', 'match_barcodes',
								  'center_pixel_x', 'center_pixel_y', 'center_pixel_z']+['%s_%d'%(mi,fr)
											    for mi,fr in zip(['max_intensity', 'second_intensity']*2*nframes, [int(d/2)
											    for d in xrange(2*nframes)])])
			    for match_barcodes in [nbr_rounds,nbr_rounds_sub]:
				    #for roi_idx in [0]:
				    for roi_idx in xrange(roi_count):
					    for barcode in patches[match_barcodes][roi_idx]:
						    for pat in patches[match_barcodes][roi_idx][barcode]:
							center_pixel = [sum([pat_pixel[0] for pat_pixel in pat['pixels']])/len(pat['pixels']),sum([pat_pixel[1] for pat_pixel in pat['pixels']])/len(pat['pixels']),sum([pat_pixel[2] for pat_pixel in pat['pixels']])/len(pat['pixels'])]
							all_intens = [pat['%s_%d'%(mi,fr)]/max_global_intensity for mi,fr in zip(['max_intensity','second_intensity']*2*nframes,[int(d/2) for d in xrange(2*nframes)])]
							tsv_writer.writerow([barcode,pat['roi'],pat['max_intensity_pixel'][0]*pixelWidth,
										      pat['max_intensity_pixel'][1]*pixelWidth,pat['max_intensity_pixel'][2]*pixelDepth,
										      pat['max_intensity'],len(pat['pixels']),pat['match_barcodes'],
										      center_pixel[0]*pixelWidth, center_pixel[1]*pixelWidth,
										      center_pixel[2]*pixelDepth]+all_intens)
print 'Finished!'
