# pylint: skip-file
import csv
import os
import getopt
from math import sqrt, exp
from collections import OrderedDict
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from copkmeans.cop_kmeans import cop_kmeans
from itertools import combinations, product
import numpy as np
import IMP
import IMP.domino
import IMP.core
import IMP.container
import sys

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

def write_overlay(directory, models, prefix='', labels=None, roi_nbr=0):

    out = 'roi,chr,x,y,z,label'
    form = ('\n%d,%s,%f,%f,%f,%s')
    for idx,model in enumerate(models):
        cur_chrom=''
        for i in xrange(len(model['x'])):
            chr = (labels[idx][i].split('p')[0]).split('q')[0] if ('p' in labels[idx][i] or 'q' in labels[idx][i]) else labels[idx][i][0] \
                    if labels else ''
            if cur_chrom != chr and chr != '':
                cur_ploid = labels[idx][i].split('-')[1]
                cur_chrom = chr + '-' + cur_ploid 
            out += form % (model['rand_init'], cur_chrom,
                           model['x'][i], model['y'][i], model['z'][i],labels[idx][i] if labels else '0')
    out_f = open('%s/%s-overlay-%d.csv' % (directory, prefix, roi_nbr), 'w')
    out_f.write(out)
    out_f.close()

def median(lst):
    sortedLst = sorted(lst)
    lstLen = len(lst)
    index = (lstLen - 1) // 2

    if (lstLen % 2):
        return sortedLst[index]
    else:
        return (sortedLst[index] + sortedLst[index + 1])/2.0
     
def create_scoring(m, ps, locis, ploid, center_clusters,
                   chr_territory, chr_join_max_dist = 3.5):
    
    restr = []
    pairs=[]
    for pl in xrange(ploid):
        for i in xrange(len(locis)-1):
            pairs.append([i+pl*len(locis), i+pl*len(locis)+1])
    pc= IMP.container.ListPairContainer(m,[(ps[p[0]], ps[p[1]]) for p in pairs],"Restrained pairs")
    score = IMP.core.DistancePairScore(IMP.core.HarmonicUpperBound(chr_join_max_dist, 1))
    pr= IMP.container.PairsRestraint(score, pc)
    pr.set_maximum_score(.01)
    restr.append(pr)
    
    for pl in xrange(ploid):
        center_cluster = center_clusters[pl]
        parts = [ps[i] for i in range(pl*len(locis),(pl+1)*len(locis))]
        sc = IMP.container.ListSingletonContainer(m,parts,"Chromosome %d"%(pl+1))
        segment = IMP.algebra.Segment3D(
                IMP.algebra.Vector3D(center_cluster[0],center_cluster[1],0),
                IMP.algebra.Vector3D(center_cluster[0],center_cluster[1],3))
        bb = IMP.algebra.get_bounding_box(segment)
        # define a restraint to keep all the model['particles'] at a distance lower or equal to the cylinder radius
        scoren = IMP.core.BoundingBox3DSingletonScore(IMP.core.HarmonicUpperBound(chr_territory,1), bb)
        prs = IMP.container.SingletonsRestraint(scoren,sc)
        prs.set_maximum_score(.01)
        restr.append(prs)
    
    return restr


def create_representation(m, locis, ploid, start_num_chr=0):
    ps = []
    for pl in xrange(ploid):
        for i, part in enumerate(locis):
            p = m.add_particle("%s-%i" % (part, start_num_chr+pl+1))
            IMP.core.XYZ.setup_particle(m, p,
                                        IMP.algebra.Vector3D([0,0,0]))
                                        #IMP.algebra.Vector3D(points[i][0]))
            ps.append(p)
    return ps


def create_discrete_states(m, ps, locis, barcode_points, center_clusters,
                           ploid, fix_locis, probs):
    
    pst = IMP.domino.ParticleStatesTable()
    assigned_points = [[] for _ in xrange(len(locis)*ploid)]   
    for i, part in enumerate(locis):
        #states = SpotStates(vs, intensities[i])
        ass_points = []
        for idx, (pb_idx, pb_i, np_point) in enumerate(barcode_points[i]):
            if pb_idx not in [fl[0] for fl in [fix_locis[i+pl*len(locis)] for pl in xrange(ploid)] if fl]:
                ass_points.append((pb_idx,np_point,probs[i][idx]))
        #points += fake_points[i]
        #ass_points += [(-1,f_point,0) for f_point in fake_points[i]]
        use_fix = False
        for pl in xrange(ploid):
            if fix_locis[i+pl*len(locis)] is not None:
                assigned_points[i+pl*len(locis)] = [fix_locis[i+pl*len(locis)]]
                use_fix = True
                #print m.get_particle(ps[i+pl*len(locis)]), fix_locis[i+pl*len(locis)][0]
            else:
                list_assigned_left = dict((i-k,fix_locis[k+pl*len(locis)][1])
                                          for k in xrange(i) if fix_locis[k+pl*len(locis)]) if i>0 else {}
                list_assigned_right = dict((k-i,fix_locis[k+pl*len(locis)][1])
                                           for k in range(i+1, len(locis)) if fix_locis[k+pl*len(locis)]) if i<len(locis)-1 else {}
                left_loci = right_loci = center_clusters[pl]
                if len(list_assigned_left) > 0:
                    left_loci = list_assigned_left[min(list_assigned_left)]
                if len(list_assigned_right) > 0:
                    right_loci = list_assigned_right[min(list_assigned_right)]
                closest_assigned = (-1,(np.array(right_loci)+np.array(left_loci))/2,0)
#                 if i+1+pl*len(locis) < len(fix_locis) and fix_locis[i+1+pl*len(locis)]:
#                     closest_assigned = fix_locis[i+1+pl*len(locis)][1]
#                 elif i-1+pl*len(locis) >= 0 and fix_locis[i-1+pl*len(locis)]:
#                     closest_assigned = fix_locis[i-1+pl*len(locis)][1]
                assigned_points[i+pl*len(locis)]= ass_points+[closest_assigned]
                
        if use_fix:
            for pl in xrange(ploid):
                points = [point[1] for point in assigned_points[i+pl*len(locis)]]
                vs = [IMP.algebra.Vector3D(point) for point in points]
                states = IMP.domino.XYZStates(vs)
                pst.set_particle_states(m.get_particle(ps[i+pl*len(locis)]), states)
                #print len(assigned_points[i+pl*len(locis)]), "states for particle ",m.get_particle(ps[i+pl*len(locis)])
                
        else:
            merged_points = [point for point in assigned_points[i] if point[0] >= 0]
            for pl in xrange(ploid):
                merged_points += [point for point in assigned_points[i+pl*len(locis)]
                                  if point and point[0] < 0]
            for pl in xrange(ploid):
                assigned_points[i+pl*len(locis)] = merged_points
            points = [point[1] for point in assigned_points[i]]
            vs = [IMP.algebra.Vector3D(point) for point in points]
            states = IMP.domino.XYZStates(vs)
            for pl in xrange(ploid):
                pst.set_particle_states(m.get_particle(ps[i+pl*len(locis)]), states)    
                #print len(assigned_points[i]), "states for particle ",m.get_particle(ps[i+pl*len(locis)])
    return pst, assigned_points


def create_sampler(m, r, pst, coord_loci):
    # create the sampler and pass it the states for each patricle
    s = IMP.domino.DominoSampler(m, pst)
    s.set_restraints(r)
    # the following lines recreate the defaults and so are optional
    filters = []
    # create a restraint cache to avoid re-evaluating restraints
    rc = IMP.domino.RestraintCache(pst)
    # add the list of restraints we want to use
    rc.add_restraints(r)
    # do not allow particles with the same ParticleStates object
    # to have the same state index
    #if len([c for c in coord_loci if c]) == 0:
    filters.append(IMP.domino.ExclusionSubsetFilterTable(pst))
    # filter states that score worse than the cutoffs in the Model
    filters.append(IMP.domino.RestraintScoreSubsetFilterTable(rc))
    filters[-1].set_log_level(IMP.SILENT)
    # try to be intelligent about enumerating the states in each subset
    states = IMP.domino.BranchAndBoundAssignmentsTable(pst, filters)
    states.set_log_level(IMP.SILENT)
    s.set_assignments_table(states)
    s.set_subset_filter_tables(filters)

    return s

def chunks(l, n):
    for i in range(0, len(l), n):
        yield l[i:i+n]

def get_constrained_clusters(best_patchs, nbr_chrs, min_cluster_separation, chr_territory,
                             max_locis_per_barcode = 4, min_loci_per_chr=2):
    
    nbr_points_barcode = []
    points = []
    for barcode in best_patchs:
        nbr_points_barcode.append(len(best_patchs[barcode][:max_locis_per_barcode]))
        points += [{'barcode':barcode,
                    'intensity': b_p['max_intensity'],
                    'point': (b_p['center_pixel_x'],
                              b_p['center_pixel_y'],
                              b_p['center_pixel_z'])}
                    for b_p in best_patchs[barcode][:max_locis_per_barcode]]
    if len(points) == 0:
        return None
    cannot_link = []
    total_points = 0
    for nbr_points in nbr_points_barcode:
        cannot_link += [comb for comb in combinations(range(total_points,total_points+nbr_points),2)]
        total_points += nbr_points
    positions_xyz = np.array([np.array(pb_i['point']) for pb_i in points]).reshape(len(points),3)
    
    detect_centers = []
    # Detect the expected clusters
    nbr_clusters = nbr_chrs
    if len(positions_xyz) >= nbr_clusters:
        clusters_ass, cluster_centers = cop_kmeans(dataset=positions_xyz, k=nbr_clusters, ml=[],cl=cannot_link)
        # Maybe we have too much spots and 
        # cannot be split
        if cluster_centers is None and len(positions_xyz) >= nbr_clusters+1:
            nbr_clusters += 1
            clusters_ass, cluster_centers = cop_kmeans(dataset=positions_xyz, k=nbr_clusters, ml=[],cl=cannot_link)
        # If not try to get at least nbr_clusters-1 
        if cluster_centers is None and len(positions_xyz) >= nbr_clusters-1:
            nbr_clusters -= 1
            clusters_ass, cluster_centers = cop_kmeans(dataset=positions_xyz, k=nbr_clusters, ml=[],cl=cannot_link)
        if cluster_centers is not None:
            for num_clust in xrange(nbr_clusters):
                barcode_clust = [pb['barcode'] for pb_i, pb in enumerate(points)
                                 if clusters_ass[pb_i] == num_clust
                                 and np.linalg.norm(cluster_centers[num_clust]-np.array(pb['point']))<chr_territory]
                intensities_clust = sum([pb['intensity'] for pb_i, pb in enumerate(points) if clusters_ass[pb_i] == num_clust])
                if len(set(barcode_clust)) >= min_loci_per_chr:
                    detect_centers.append({'intensity': intensities_clust/len(set(barcode_clust)),
                                           'nbr_barcodes': len(set(barcode_clust)), 
                                           'center':np.array(cluster_centers[num_clust])})
                
    detect_centers = sorted(detect_centers,key=lambda k: (k['nbr_barcodes'],k['intensity']),reverse=True)
    
    optim_centers = []
    for clust_cent in detect_centers:
        same_position = False
        for op_clust_cent in xrange(len(optim_centers)):
            if np.linalg.norm(clust_cent['center']-optim_centers[op_clust_cent]) < min_cluster_separation:
                optim_centers[op_clust_cent] = (clust_cent['center']+optim_centers[op_clust_cent])/2
                same_position = True
                break
        if not same_position:
            optim_centers.append(clust_cent['center'])
        
    #if len(optim_centers) == 0:
    #    optim_centers = [np.mean(positions_xyz, axis=0)]
    return optim_centers

if __name__ == "__main__":
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hi:r:c:b:n:s:t:o:",["puncta_file=",
                                                                   "roi=",
                                                                   "chrs=",
                                                                   "barcode_file=",
                                                                   "rounds=",
                                                                   "shared_barcodes=",
                                                                   "threshold=",
                                                                   "outdir="])
    except getopt.GetoptError:
        print '''domino_oligo_cluster_chrX.py -i --puncta_file -r --roi 
            -c --chrs -b --barcode_file -n --rounds -s --shared_barcodes
            -t --threshold -o --outdir'''
        sys.exit(2)
        
    chrs = roi = selected_roi = barcode_file = nbr_rounds  = outdir = None
    shared_barcodes = []
    
    for opt, arg in opts:
        if opt == '-h':
            print '''domino_oligo_cluster.py -i --puncta_file -r --roi 
            -c --chrs -p --plex -o --outdir'''
            sys.exit()
        elif opt in ("-i", "--puncta_file"):
            puncta_file = arg
        elif opt in ("-r", "--roi"):
            roi = arg.split(' ')
        elif opt in ("-c", "--chrs"):
            chrs = arg.split(' ')
        elif opt in ("-b", "--barcode_file"):
            barcode_file = arg
        elif opt in ("-n", "--rounds"):
            nbr_rounds = int(arg)
        elif opt in ("-s", "--shared_barcodes"):
            shared_barcodes = [int(a) for a in arg.split(',')]
        elif opt in ("-t", "--threshold"):
            threshold = [str(a) for a in arg.split(',')]
        elif opt in ("-o", "--outdir"):
            outdir = arg

    puncta_filename = os.path.splitext(os.path.basename(puncta_file))[0]
    
    probes = load_probes(barcode_file)
    barcodes = [int(probe['Readout']) for probe in probes]
    locis = [probe['Probe'] for probe in probes]
    chrs = chrs or list(set([(probe['Chr.'].split('p')[0]).split('q')[0] for probe in probes]))
    
    print puncta_file, chrs, roi, barcode_file, outdir, locis, chrs
    
    nbr_rounds = nbr_rounds or 5
    
    nbr_barcode_subsample = 4
    min_nbr_pixels = 1
    #min_nbr_pixels_best = 5
    min_nbr_pixels_best = 2
    #min_loci_per_chr = 6
    min_loci_per_chr = 3
    
    min_cluster_separation = 1.
    max_best_locis = 6
    threshold_mitotic = 2.5
    chr_join_max_dist_46plex = 2.5
    chr_join_max_dist_23plex = 3
    #chr_territory = 4.5
    chr_territory = 3.5
    pixelWidth = 0.267
    pixelDepth = 0.3
    ImageWidth = 512
    ImageHeight = 512
    max_points_domino = 4
    barcode_step = 20
    thresholds = {('JEB','5','MS','True'):(62000./(2**16-1),61000./(2**16-1),47000),
                  ('JEB','5','MS','False'):(62000./(2**16-1),61000./(2**16-1),0),
                  ('JEB','5','MSBS','True'):(64000./(2**16-1),63500./(2**16-1),47000),
                  ('JEB','5','MSBS','False'):(64000./(2**16-1),63500./(2**16-1),0),
                  ('SOLiD','4','MSBS','False'):(64000./(2**16-1),63500./(2**16-1),0),
                  ('SOLiD','4','MSBS','True'):(64000./(2**16-1),63500./(2**16-1),47000),
                  ('SOLiD','5','MSBS','False'):(64000./(2**16-1),63500./(2**16-1),0),
                  ('SOLiD','4','MS','False'):(59000./(2**16-1),58000./(2**16-1),0),
                  ('IMR90','5','MSBS','False'):(61000./(2**16-1),60000./(2**16-1),0)}
    key_location = 'max_intensity'
    
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    
    max_int_label = 'max_intensity'
    if threshold[0] == 'SOLiD':
        max_int_label = 'avg_intensity3'
    sel_threshold = thresholds[tuple(threshold)]
    
    patches = OrderedDict((barc,[]) for barc in barcodes)
    totos = {}
    file_rois = []
    if sel_threshold[2] > 0:
        with open(puncta_file.replace('_punctas','_toto_intensities'), 'r') as csvfile:
            csvreader = csv.reader(csvfile, delimiter='\t')
            header = next(csvreader, None)
            for row in csvreader:
                totos[(int(row[0]),int(row[1]),int(row[2]),int(row[3]))] = float(row[4])
    #with open(input_folder+'punctas.csv', 'r') as csvfile:
    #with open(input_folder+'cell 4/merged/OFQv41_w4_Chr235_newMSBS_setA_001_all_punctas.tsv', 'r') as csvfile:
    with open(puncta_file, 'r') as csvfile:
    #with open(input_folder+'raw_nd2/merged/OFQv41_w4_Chr235_newMSBS_setA_001_punctas.tsv', 'r') as csvfile:
        csvreader = csv.reader(csvfile, delimiter='\t')
        header = next(csvreader, None)
        for row in csvreader:
            spot = {}
            for i,h in enumerate(header):
                if '.' in row[i]:
                    spot[h] = float(row[i])
                else:
                    spot[h] = int(row[i])
            pixel_x, pixel_y, pixel_z = int(round(spot['max_intensity_pixel_x']/pixelWidth)), \
                int(round(spot['max_intensity_pixel_y']/pixelWidth)),int(round(spot['max_intensity_pixel_z']/pixelDepth))
            spot['avg_intensity3'] = (spot['max_intensity_0']+spot['max_intensity_1']+spot['max_intensity_2'])/3.
            spot['toto_intensity'] = 0
            if (spot['roi'],pixel_x, pixel_y, pixel_z) in totos:
                spot['toto_intensity'] = totos[(spot['roi'],pixel_x, pixel_y, pixel_z)]
            if spot['roi'] not in file_rois:
                file_rois.append(spot['roi'])
            if spot['barcode'] in patches:
                #Filter coincident
                #if spot['match_barcodes'] <  nbr_rounds:
                spot_loc = np.array([spot[key_location+'_pixel_x'],
                                     spot[key_location+'_pixel_y']])
                same_spots = [pb for pb in patches[spot['barcode']] if np.linalg.norm(np.array([pb[key_location+'_pixel_x'],
                                                                                                pb[key_location+'_pixel_y']]) - spot_loc) < 1]
                if len(same_spots) > 0:
                    continue
                patches[spot['barcode']].append(spot)
            
    all_intensities = []
    all_areas = []
    all_qualities = []
    for barcode in barcodes:
        patches[barcode] = sorted(patches[barcode],key=lambda k: k[max_int_label], reverse=True)
        all_intensities += [fp[max_int_label] for fp in patches[barcode] if fp['match_barcodes']==nbr_rounds]
        all_areas += [fp['nbr_pixels'] for fp in patches[barcode] if fp['match_barcodes']==nbr_rounds]
    #print all_qualities
    median_intensity = max(all_intensities)*0.75
    min_intensity = median_intensity
    min_intensity_single_pixel = median_intensity*1.1
    #min_intensity_single_pixel = median_intensity*1.5
    print 'Using filter by minimum intensity:', min_intensity, min_intensity_single_pixel
    
    all_labels = []
    models = []
    #roi = [1]
    if not roi:
        selected_roi = file_rois
    else:
        selected_roi = [int(r) for r in roi]
    all_nbr_detected_loci = 0
    output_file = open('%s/%s_out.tsv' % (outdir, puncta_filename), mode='w')
    tsv_writer = csv.writer(output_file, delimiter='\t', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    tsv_writer.writerow(['roi', 'chr', 'loci', 'pos_x', 'pos_y', 'pos_z', 'score',
                         'nbr_pixels', 'max_intensity', 'match_barcodes'])
    all_best_patches = {}
    #for roi_idx in [0,1,2,3,4,5,6]:
    #chrs = ['3']
    for roi_idx in selected_roi:
        print 'Processing roi', roi_idx+1
        labels = []
        nbr_detected_loci = 0
        total_possible_locis = 0
                        
        chrom_centers = []
        detected_chrs = []
        avg_barcodes = 0
        for chrom in chrs:
            used_chr_territory = chr_territory
            chrom_barcodes = [bar for idx_b,bar in enumerate(barcodes) if locis[idx_b][:len(chrom)] == chrom]
            best_patchs = {}
            for barcode in chrom_barcodes:
                if barcode in shared_barcodes:
                    best_patchs[barcode] = []
                    continue
                    
                best_patchs[barcode] = [fp for fp in patches[barcode] if fp['roi'] == roi_idx
                                    and fp['match_barcodes']==nbr_rounds
                                    and fp[max_int_label]>sel_threshold[1]
                                    and (fp['toto_intensity']>sel_threshold[2] or sel_threshold[2] == 0)
                                    and fp['nbr_pixels']>min_nbr_pixels_best]
                best_patchs[barcode] = sorted(best_patchs[barcode],key=lambda k: k['max_intensity'],reverse=True)
                if len(best_patchs[barcode]) > 0:
                    avg_barcodes = (avg_barcodes + len(best_patchs[barcode]))
            
            all_best_patches.update(best_patchs)
        
        # decide with if this cell is in M phase
        ploidy = 1
        diploid_chroms = [chrom for chrom in chrs if chrom != 'X']
        if len(diploid_chroms) > 0:
            if avg_barcodes > threshold_mitotic:
                ploidy = 2
        else:
            if avg_barcodes > threshold_mitotic/2.:
                ploidy = 2

        for chrom_idx, chrom in enumerate(chrs):
            chrom_barcodes = [bar for idx_b,bar in enumerate(barcodes) if locis[idx_b][:len(chrom)] == chrom]
            max_nbr_chrs = 2*ploidy
            if chrom == 'X':
                max_nbr_chrs = 1*ploidy
            nbr_max_best_locis = max_nbr_chrs
            
            best_patchs = {barcode: all_best_patches[barcode] for barcode in chrom_barcodes if barcode not in shared_barcodes}
            c_cluster = get_constrained_clusters(best_patchs, max_nbr_chrs, min_cluster_separation, 
                                                 chr_territory, max_locis_per_barcode=nbr_max_best_locis,
                                                 min_loci_per_chr=min_loci_per_chr)
            if c_cluster is not None:
                detected_chrs.append(chrom)
                chrom_centers.append(c_cluster)
                print 'Centers for chromosome ', chrom, c_cluster
            #c_cluster = []    
            #for _ in xrange(nbr_chrs):
            #    c_cluster += get_centers(chrom_barcodes, nbr_chrs, chr_territory, min_cluster_separation,
            #                            nbr_max_best_locis, min_loci_per_chr=min_loci_per_chr,
            #                            exclude_centers=[(c_[0],c_[1]) for c_ in c_cluster])
                            
        if len(chrom_centers) == 0:
            continue

        # just skip mitotic for the 36plex for the time being
        if ploidy > 1 and len(chrom_centers) > 2:
            print("Skipping cell in m-phase")
            continue

        models.append({'rand_init':roi_idx,'x':[],'y':[],'z':[]})
        model = models[-1]
        for chrom_idx, chrom in enumerate(detected_chrs):
            center_clusters=chrom_centers[chrom_idx]
            if len(center_clusters) == 0:
                continue
            used_chr_territory = chr_territory
            used_max_points_domino = max_points_domino
            if len(center_clusters) > 2:
                used_max_points_domino = int(max_points_domino/2)
            chrom_barcodes = [bar for idx_b,bar in enumerate(barcodes) if locis[idx_b][:len(chrom)] == chrom]

            chr_join_max_dist = chr_join_max_dist_23plex
            if len(chrom_barcodes) > 23:
                chr_join_max_dist = chr_join_max_dist_46plex
                
            used_nbr_chrs = len(center_clusters)
            #chrom_intensities = [fp['max_intensity'] for barcode in chrom_barcodes for fp in patches[barcode]]
            #chrom_nbr_pixels = [fp['nbr_pixels'] for barcode in chrom_barcodes for fp in patches[barcode]]            
            #chrom_median_intensity = median(chrom_intensities)
            #chrom_median_intensity = sum(chrom_intensities)/len(chrom_intensities)
            #chrom_max_intensity = max(chrom_intensities)
            #intensity_thresholds = list(np.linspace(start = chrom_max_intensity, stop = min_intensity, num = intensity_steps))
            #nbr_pixels_thresholds = list(np.linspace(start = max(chrom_nbr_pixels), stop = 1, num = intensity_steps))
            coord_loci = None
            for barcode_category in ['strong positives','full barcodes','subsampling']:
                roi_intensities = []
                used_barcode_step = barcode_step
                min_intensity_threshold = sel_threshold[1]
                min_toto_threshold = sel_threshold[2]
                match_barcodes_threshold = nbr_rounds
                if barcode_category == 'full_barcodes':
                    min_intensity_threshold = min_intensity
                    min_toto_threshold = 0
                elif barcode_category == 'subsampling':
                    match_barcodes_threshold = nbr_barcode_subsample
                    min_intensity_threshold = min_intensity
                    min_toto_threshold = 0
                
                for barcode in chrom_barcodes:
                    roi_intensities += [{'id': fp_id, 'barcode': barcode, 'max_intensity': fp['max_intensity']}
                                        for fp_id, fp in enumerate(patches[barcode]) if fp['roi'] == roi_idx
                                        and fp['match_barcodes']>=match_barcodes_threshold
                                        and (fp['toto_intensity']>min_toto_threshold or min_toto_threshold == 0)
                                        and (fp['max_intensity']>min_intensity_threshold)
                                        and (fp['nbr_pixels'] > 2 or fp['max_intensity']>min_intensity_single_pixel)
                                        ]
                if len(roi_intensities) == 0:
                    continue
                roi_intensities = sorted(roi_intensities,
                                          key=lambda k: k['max_intensity'],reverse=True)
                used_barcodes = 0
                while used_barcodes < len(roi_intensities):
                #for chrom_intensity_threshold, nbr_pixels_threshold in zip(intensity_thresholds,nbr_pixels_thresholds):
                    #print 'Using matches %d and threshold %f' % (match_barcodes_threshold,
                    #                                                 chrom_intensity_threshold)
                    patches_ids = {}
                    for barcode in chrom_barcodes:
                        patches_ids[barcode] = [pat['id']
                                                for pat in roi_intensities[used_barcodes:used_barcodes+used_barcode_step]
                                                if pat['barcode'] == barcode]
                        patches_ids[barcode] = sorted(patches_ids[barcode],
                                                      key=lambda k: (patches[barcode][k]['match_barcodes'],
                                                                     patches[barcode][k]['nbr_pixels']),
                                                      reverse=True)[:used_max_points_domino]
                    
                    used_barcodes += min(used_barcode_step,len(roi_intensities)) 
                    if sum([len(patches_ids[barcode]) for barcode in chrom_barcodes]) == 0:
                        continue
                        
                    chrom_locis = [loci for loci in locis if loci[:len(chrom)] == chrom]
                    points_barcodes = []
                    intensities = []
                    points = []
                    if not coord_loci:
                        coord_loci = [None]*used_nbr_chrs*len(chrom_locis)
                    
                    IMP.set_log_level(IMP.SILENT)
                    m = IMP.Model()
                    # don't print information during Model.evaluate
                    m.set_log_level(IMP.SILENT)
                    
                    all_points_barcode = None
                    #fake_points = [[] for _ in xrange(len(chrom_barcodes))]
                    used_spots_ids = dict((barcode,[]) for barcode in chrom_barcodes)
                    for clust_num,clus in enumerate(center_clusters):
                        points_barcodes = []
                        for barcode in chrom_barcodes:
                            points_barcode = []
                            used_spots_id = []
                            for b_id in patches_ids[barcode]:
                                if b_id in used_spots_ids[barcode]:
                                    continue
                                pb_i = patches[barcode][b_id]
                                key_location = 'max_intensity'
                                if np.linalg.norm(np.array([pb_i[key_location+'_pixel_x'],
                                                            pb_i[key_location+'_pixel_y'],
                                                            pb_i[key_location+'_pixel_z']]) - clus) < used_chr_territory:
                                    points_barcode.append((b_id, pb_i, np.array([pb_i[key_location+'_pixel_x'],
                                                                                 pb_i[key_location+'_pixel_y'],
                                                                                 pb_i[key_location+'_pixel_z']])))
                                    used_spots_ids[barcode].append(b_id)
                            points_barcodes.append(points_barcode) 
                        
                        #for id_p,barcode in enumerate(chrom_barcodes):
                        #    fake_points[id_p] += [clus]
                            
                        
                        if all_points_barcode:
                            all_points_barcode = [all_point_barcode+point_barcode
                                                  for all_point_barcode,point_barcode in zip(all_points_barcode,points_barcodes)]
                        else:
                            all_points_barcode = points_barcodes
                    
                    #for id_p,barcode in enumerate(chrom_barcodes):         
                    #    all_points[id_p] = all_points[id_p]+fake_points[barcode]
                    #intensities = [[pb_i[1]['max_intensity'] for pb_i in points_barcode] 
                    #               for points_barcode in all_points_barcode]
                    #areas = [[float(pb_i[1]['nbr_pixels'])/max(all_areas)+float(pb_i[1]['max_intensity'])/max(all_intensities)
                    #          for pb_i in points_barcode] 
                    #               for points_barcode in all_points_barcode]
                    probs = [[
                              (pb_i[1]['nbr_pixels']-min(all_areas))/max(all_areas) + \
                              (float(pb_i[1]['max_intensity'])-min(all_intensities))/max(all_intensities)
                              for pb_i in points_barcode] 
                                   for points_barcode in all_points_barcode]
                    #probs = [[10*(pb_i[1]['spot_quality']-min(all_qualities))/max(all_qualities) + \
                    #          pb_i[1]['nbr_pixels'] + \
                    #          5*(float(pb_i[1]['max_intensity'])-min(all_intensities))/max(all_intensities)
                    #          for pb_i in points_barcode] 
                    #               for points_barcode in all_points_barcode]
                    #print "creating representation"
                    ps = create_representation(m, chrom_locis, used_nbr_chrs)
                    #print "creating discrete states"
                    pst, assigned_points = create_discrete_states(m, ps, chrom_locis, all_points_barcode, center_clusters,
                                                                  used_nbr_chrs, coord_loci, probs)
                    #print "creating score function"
                    rs = create_scoring(m, ps, chrom_locis, used_nbr_chrs, center_clusters, used_chr_territory,
                                        chr_join_max_dist=chr_join_max_dist)
                    #print "creating sampler"
                    s = create_sampler(m, rs, pst, coord_loci)
                
                    #print "sampling"
                    
                    subset = pst.get_subset()
                    #print "ordering subset"
                    ps_subset = sorted(range(len(subset)), key=lambda k: subset[k].get_index())
                    print "getting solutions"
                    solutions = s.get_sample_assignments(subset)
                    
                    if len(solutions) == 0:
                        print("There are no solutions to the problem")
                    else:
                        print " There are", len(solutions), "possible solutions"
                        best_solution = []
                        for assignment_chunk in chunks(solutions,100000):
                            ass_prob = []
                            for assignment in assignment_chunk:
                                ass_points = []
                                for pl_n in xrange(used_nbr_chrs):
                                    for idx in xrange(len((chrom_locis))):
                                        ass_points.append(assigned_points[idx+pl_n*len(chrom_locis)][assignment[ps_subset[idx+pl_n*len(chrom_locis)]]])
                                #ass_points = [assigned_points[i+pl*len(chrom_locis)][assignment[ps_subset[i+pl*len(chrom_locis)]]] 
                                #              for pl in xrange(nbr_chrs) for i in xrange(len((chrom_locis)))]
                                #if len(assigned_points[i+pl*len(chrom_locis)]) > assignment[ps_subset[i+pl*len(chrom_locis)]]]
                                total_length = 0
                                for pl_n in xrange(used_nbr_chrs):
                                    real_ass_points = [ass_points[k][1] for k in range(pl_n*len(chrom_locis),(pl_n+1)*len(chrom_locis)) if ass_points[k][0]>=0]
                                    for p in xrange(len(real_ass_points)-1):
                                        total_length += np.linalg.norm(real_ass_points[p+1]-real_ass_points[p])
                                
                                #total_prob = sum([intensities[i][assignment[ps_subset[i+pl*len(chrom_locis)]]]
                                #                  if len(intensities[i]) > assignment[ps_subset[i+pl*len(chrom_locis)]] else 0 
                                #                  for i in xrange(len((chrom_locis))) for pl in xrange(nbr_chrs)])
                                total_prob = sum([assigned_points[i+pl*len(chrom_locis)][assignment[ps_subset[i+pl*len(chrom_locis)]]][2]
                                                  #if len(assigned_points[i+pl*len(chrom_locis)]) > assignment[ps_subset[i+pl*len(chrom_locis)]] else 0
                                                  for i in xrange(len((chrom_locis))) for pl in xrange(used_nbr_chrs)])
                                ass_prob.append({'assignment':assignment,'prob':total_prob, 'length': total_length})
                            ass_prob = sorted(ass_prob, key=lambda k: k['prob'],reverse=True)
                            most_probable = ass_prob[0]['prob']
                            ass_prob = [ap for ap in ass_prob if ap['prob']==most_probable]
                            ass_prob = sorted(ass_prob, key=lambda k: k['length'])
                            best_solution.append(ass_prob[0])
                        
                        ass_prob = best_solution
                        ass_prob = sorted(ass_prob, key=lambda k: k['prob'],reverse=True)
                        most_probable = ass_prob[0]['prob']
                        ass_prob = [ap for ap in ass_prob if ap['prob']==most_probable]
                        ass_prob = sorted(ass_prob, key=lambda k: k['length'])
                        best_solution = ass_prob[0]['assignment']

                        print "=================> BEST SOLUTION <=============="
                        print "prob", most_probable
                        for pl_n in xrange(used_nbr_chrs):
                            for id_p, part in enumerate(chrom_locis):
                                idx_subset = ps_subset[id_p+pl_n*len(chrom_locis)]
                                if coord_loci[id_p+pl_n*len(chrom_locis)] is not None \
                                                or assigned_points[id_p+pl_n*len(chrom_locis)][best_solution[idx_subset]][0] < 0:
                                    continue
                                #coord_loci[id_p+pl_n*len(chrom_locis)] = (all_points[id_p][best_solution[idx_subset]],probs[id_p][best_solution[idx_subset]])
                                coord_loci[id_p+pl_n*len(chrom_locis)] = (assigned_points[id_p+pl_n*len(chrom_locis)][best_solution[idx_subset]][0],
                                                                          assigned_points[id_p+pl_n*len(chrom_locis)][best_solution[idx_subset]][1],
                                                                          assigned_points[id_p+pl_n*len(chrom_locis)][best_solution[idx_subset]][2])
                                label_loci = subset[idx_subset]
                                print label_loci, coord_loci[id_p+pl_n*len(chrom_locis)]
                    if len([coord_loc for coord_loc in coord_loci if coord_loc is None]) == 0:
                        break
                                    
    #                         for pl_n in xrange(nbr_chrs):
    #                             for id_p, part in enumerate(chrom_locis):
    #                                 idx_subset = ps_subset[id_p+pl_n*len(chrom_locis)]
    #                                 all_points[id_p][(best_solution[idx_subset])] = None
    #                                 if len(probs[id_p]) <= best_solution[idx_subset]:
    #                                     continue
    #                                 all_points_barcode[id_p][(best_solution[idx_subset])] = None
    #                         all_points = [[all_point for all_point in all_points_bar if all_point is not None]
    #                                       for all_points_bar in all_points]
    #                         all_points_barcode = [[all_point for all_point in all_points_bar if all_point is not None]
    #                                       for all_points_bar in all_points_barcode]
                    print '%d/%d' % (used_barcodes, len(roi_intensities))
            #Esto dentro del bucle
            for pl_n in xrange(used_nbr_chrs):
                total_possible_locis += len(chrom_locis)
                for id_p, part in enumerate(chrom_locis):
                    idx_subset = ps_subset[id_p+pl_n*len(chrom_locis)]
                    if coord_loci[id_p+pl_n*len(chrom_locis)] is None:
                        continue
                    #coord_loci[id_p+pl_n*len(chrom_locis)] = all_points[id_p][best_solution[idx_subset]]
                    label_loci = subset[idx_subset] 
                    model['x'].append(coord_loci[id_p+pl_n*len(chrom_locis)][1][0])
                    model['y'].append(coord_loci[id_p+pl_n*len(chrom_locis)][1][1])
                    model['z'].append(coord_loci[id_p+pl_n*len(chrom_locis)][1][2])
                    labels.append(str(label_loci).replace('"',''))
                    nbr_detected_loci += 1
                    
                    #print label_loci, coord_loci[id_p+pl_n*len(chrom_locis)][1]
                    lab_roi_idx = str(roi_idx) if pl_n == 0 else str(roi_idx)+'-'+str(pl_n+1)
                    row_loci = [lab_roi_idx, chrom, str(label_loci).replace('"','')]+list(coord_loci[id_p+pl_n*len(chrom_locis)][1]) + \
                        [coord_loci[id_p+pl_n*len(chrom_locis)][2],
                         patches[chrom_barcodes[id_p]][coord_loci[id_p+pl_n*len(chrom_locis)][0]]['nbr_pixels'],
                         patches[chrom_barcodes[id_p]][coord_loci[id_p+pl_n*len(chrom_locis)][0]]['max_intensity'],
                         patches[chrom_barcodes[id_p]][coord_loci[id_p+pl_n*len(chrom_locis)][0]]['match_barcodes']]
                    tsv_writer.writerow(row_loci)
                    print row_loci
        all_labels.append(labels)
        if total_possible_locis > 0:
            print 'Percentage of detected loci for roi %d, %f' % (roi_idx, 100*(float(nbr_detected_loci)/float(total_possible_locis)))
            all_nbr_detected_loci += (float(nbr_detected_loci)/float(total_possible_locis))
    output_file.close()
    print 'Percentage of detected loci, %f' % (100*(float(all_nbr_detected_loci)/len(selected_roi)))                
    write_overlay(directory=outdir, models=models, prefix=puncta_filename, labels=all_labels, roi_nbr=roi_idx)
