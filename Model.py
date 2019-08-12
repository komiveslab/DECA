'''
2019-02-11
Two additional data formats added to import
Program remembers correction data
Fully Deuterated controls can be used for back exchange correction
Peptide RT prediction added
Pymol script export can now be Relative Uptake, Relative SD, Uptake Diff, Fractional Uptake, RSD, Fractional Uptake Difference

Early 2019
Updated to be compatible with new file formats

2019-04-19
Pymol Scripts updated to color residues without coverage black
Scripts also are specific to chain now, allowing for multiple spectrum IDs.
Import DnX files
Generate 2D and 3D spectral maps
Generate State Data from DnX

2019-04-29
Implemented Statistical Significance Evaluation

2019-07-16
'''

import sqlite3
import csv
from numpy import mean, std, sqrt, abs, array, percentile, max, min, zeros, sum, nan, ones, triu_indices
from copy import copy, deepcopy
from itertools import islice, combinations
from os.path import split, splitext
from operator import itemgetter
from pyteomics import achrom
from statsmodels.stats.weightstats import DescrStatsW, ttest_ind
from scipy.stats import levene,f_oneway, ttest_ind_from_stats
from statsmodels.stats.multicomp import MultiComparison
from statsmodels.stats.libqsturng import psturng, qsturng
from scipy.stats import f

class DECA:
    '''
    Model (Data Management) for DECA Program
    Stores data from imported files
    Corrects, Modifies, and Exports data
    '''

    # Initialize object, trigger file detection and read if not empty
    def __init__(self, name='', empty=False, que=None):
        self.master_csv = []
        self.seq_nest = {}
        self.proteins = []
        self.states = []
        self.exposures = []
        self.residues = {}
        self.sequences = {}
        self.data_nest = {}
        self.cor_factor = {}
        self.rettimes = []
        self.corrected = 'No'
        self.recombined = False
        self.merge_and_corrected = False
        self.replicates = {}
        self.stats_nest = None
        self.pepid = 0
        self.que = que
        if not empty:
            self.filename = name
            self.detectType(name)
            self.organizeData()
            self.cleanFiles()
            if self.filetype == 'DNX':
                self.assignIDs_DNX()
            else:
                self.assignIDs_CSV()
            if self.que is not None:
                self.que.put(-1)

    def dictify(self, order):
        newdict = {}
        for row in self.master_csv:
            here = newdict
            for item in order[:-1]:
                if row[item] not in here:
                    here[row[item]] = {}
                here = here[row[item]]
            if row[order[-1]] not in here:
                here[row[order[-1]]] = []
            here[row[order[-1]]].append(row)
        return newdict

    def sort(self,order):
        newlist = self.master_csv
        for item in order[::-1]:
            newlist = sorted(newlist,key=itemgetter(item))
        return newlist

    # Detects the file type based on the extension and the header of the file
    def detectType(self, infile):
        print('Detecting File')
        if str(splitext(split(infile)[1])[1]).lower() == '.dnx':
            self.aquireIons()
            self.filetype = 'DNX'
            self.readFile(infile, 'DNX')
            if self.que is not None:
                self.que.put(100)
        else:
            with open(infile, 'r') as csvfile:
                csvfile.seek(0)
                has_header = csv.Sniffer().has_header(''.join([line.replace(" ", "") for line in csvfile]))
                csvfile.seek(0)
                header = list(islice(csvfile, 0, 1))[0].split(',')
                if has_header and ('Start' in header or 'start' in header):
                    n = 0
                else:
                    for n in range(1, 100):
                        csvfile.seek(0)
                        print(n)
                        f = [line.replace(" ", "") for line in islice(csvfile, n, 100)]
                        print(f[0])
                        if f[0] == '' or f[1] == '' or f[0].count(',')<=1:
                            continue
                        r = ''.join(f)
                        has_header = csv.Sniffer().has_header(r)
                        if has_header:
                            csvfile.seek(0)
                            header = list(islice(csvfile, n, n + 1))[0].split(',')
                            break
            csvfile.close()
            if n < 99:
                header_position = n
                if all(item in header for item in ['Exposure', 'Uptake']):
                    print('File is DynamX Style')
                    self.filetype = 'CSV'
                    self.readFile(infile, 'Waters')
                elif all(item in header for item in ['timepoint', 'percentd_replicate']):
                    print('File is HDXWorkbench Style')
                    self.readFile(infile, 'HDX', header_position)
                    self.filetype = 'CSV'
                elif any(item in header for item in ['Uptake_0','Uptake_0s','Uptake_0m','Uptake_0h']):
                    print('File is Alternate Style')
                    self.filetype = 'CSV'
                    self.readFile(infile, 'ALT')
                else:
                    print('File is Unrecognized')
            else:
                print('Header not found')

    # Reads data into master_csv
    def readFile(self, infile, ftype, header_position=None):
        if ftype == 'Waters':
            exp_factors_set = set()
            self.pep_corrected = False
            row_count = sum(1 for row in open(infile))
            with open(infile, 'r') as csvfile:
                reader = csv.DictReader(csvfile)
                count = 1
                print(row_count)
                for row in reader:
                    if self.que is not None:
                        self.que.put(count/(row_count-1) * 100)
                    row['size'] = int(row['End']) - int(row['Start']) + 1
                    row['Start'] = int(row['Start'])
                    row['End'] = int(row['End'])
                    row['Protein'] = str(row['Protein'])
                    if 'Fragment' in row.keys():
                        if row['Fragment'] != '':
                            self.recombined = True
                    else:
                        row['Fragment'] = str('')
                    up_count = 0
                    for char in list(row['Sequence'])[1:]:
                        if char != "P":
                            up_count += 1
                    row['MaxUptake'] = int(up_count)
                    count = count + 1
                    self.master_csv.append(row)
                    if 'Uptake_corr' in reader.fieldnames:
                        if 'CorrectionFactor' in reader.fieldnames:
                            row['PeptideCorrection'] = 1
                            row['ExposureCorrection'] = row['CorrectionFactor']
                        elif 'ExposureCorrection' in reader.fieldnames:
                            if float(row['PeptideCorrection']) != 1:
                                self.pep_corrected = True
                        else:
                            row['ExposureCorrection'] = round(float(row['Uptake'])/float(row['Uptake_corr']),3)
                        exp_factors_set.add((row['Exposure'], row['ExposureCorrection']))
                    self.rettimes.append((row['Sequence'],row['RT']))
            if 'Uptake_corr' in list(self.master_csv[0].keys()):
                self.corrected = 'Yes'
                if any(item != 1 for item in [i[1] for i in exp_factors_set]):
                    self.exp_corrected = True
            exp_factors_list = [item[0] for item in list(exp_factors_set)]
            if len(set([x for x in exp_factors_list if exp_factors_list.count(x) > 1]))>0:
                self.merge_and_corrected = True
            for i in sorted(list(exp_factors_set)):
                self.cor_factor[float(i[0])] = float(i[1])
            rettimes2 = [i for i in self.rettimes if i[1] != '']
            self.seq_to_rt = achrom.get_RCs([str(i[0]) for i in rettimes2], [float(i[1]) for i in rettimes2],rcond=None)
            print('read finished')

        elif ftype == 'HDX':
            exp_factors_set = set()
            self.pep_corrected = False
            with open(infile, 'rU') as csvfile:
                csvfile.seek(0)
                reader = csv.DictReader(islice(csvfile, header_position, None))
                for row in reader:
                    if row['timepoint'] == 'NA':
                        continue
                    row['Sequence'] = str(row['peptide'])
                    row['Start'] = int(row['start'])
                    row['End'] = int(row['end'])
                    row['size'] = row['End'] - row['Start'] + 1
                    row['Protein'] = str(row['project'])
                    row['Exposure'] = round(float(row['timepoint'][:-1]) / 60, 2)
                    row['State'] = str(row['sample'])
                    row['RT'] = [(float(row['rt_end_replicate']) + float(row['rt_start_replicate'])) / 2]
                    row['RT SD'] = float()
                    row['Center'] = [float(row['centroid']) * float(row['charge'])]
                    row['Center SD'] = float()
                    row['Fragment'] = str('')
                    row['Modification'] = str(row['features'])
                    up_count = 0
                    for char in list(row['Sequence'])[1:]:
                        if char != "P":
                            up_count += 1
                    row['MaxUptake'] = int(up_count)
                    row['Centroid_Uptake'] = [(float(row['centroid']) - float(row['mono_mz0'])) * float(row['charge'])]
                    row['Centroid_Uptake_SD'] = float()
                    row['Envelope_Uptake'] = [float(row['MaxUptake']) * float(row['percentd_replicate'])]
                    row['Envelope_Uptake_SD'] = float()
                    row['Uptake'] = float()
                    row['Uptake SD'] = float()
                    row['MHP'] = float(row['monoisotopic']) + 1.007

                    matching_row = [i for i in self.master_csv if
                                    (row['Start'] == i['Start']) and (row['End'] == i['End'])
                                    and (row['Protein'] == i['Protein']) and (row['State'] == i['State']) and
                                    (row['MHP'] == i['MHP']) and (row['Fragment'] == i['Fragment']) and
                                    (row['Exposure'] == i['Exposure'])]
                    if len(matching_row) > 1:
                        print('Error: read_hdxworkbench, duplicate entries')
                    elif len(matching_row) == 1:
                        matching_row[0]['Center'].extend(row['Center'])
                        matching_row[0]['RT'].extend(row['RT'])
                        matching_row[0]['Centroid_Uptake'].extend(row['Centroid_Uptake'])
                        matching_row[0]['Envelope_Uptake'].extend(row['Envelope_Uptake'])
                    else:
                        self.master_csv.append(row)

                    if 'Uptake_corr' in reader.fieldnames:
                        exp_factors_set.add((row['Exposure'], row['ExposureCorrection']))
                        if float(row['PeptideCorrection']) != 1:
                            self.pep_corrected = True

                if 'Uptake_corr' in reader.fieldnames:
                    self.corrected = 'Yes'
                    if any(item != 1 for item in [i[1] for i in exp_factors_set]):
                        self.exp_corrected = True

            for row in self.master_csv:
                row.update({'Center SD':std(row['Center']),
                            'Center':mean(row['Center']),
                            'RT SD':std(row['RT']),
                            'RT':mean(row['RT']),
                            'Centroid_Uptake_SD':std(row['Centroid_Uptake']),
                            'Centroid_Uptake':mean(row['Centroid_Uptake']),
                            'Envelope_Uptake_SD':std(row['Envelope_Uptake']),
                            'Envelope_Uptake':mean(row['Envelope_Uptake']),
                            'Uptake SD': round(row['Centroid_Uptake_SD'], 2),
                            'Uptake':round(row['Centroid_Uptake'], 2)})

            for i in sorted(list(exp_factors_set)):
                self.cor_factor[float(i[0])] = float(i[1])

        elif ftype == 'DNX':
            states = list(set([i['State'] for i in self.rawfiles]))
            exposures = list(set([float(i['Exposure']) for i in self.rawfiles]))
            master_csv = []
            n = 0
            exp_factors_set = set()
            zeroes = {}
            for peptide in self.peptides:
                if self.que is not None and round((n/len(self.peptides)*75)+25,0)%5==0:
                    self.que.put(round((n/len(self.peptides)*75)+25,0))
                with sqlite3.connect(self.filename) as conn:
                    conn.row_factory = sqlite3.Row
                    c = conn.cursor()
                    c.execute(
                        '''SELECT c.ID as PepID, d.name as Protein, f.name as State, c.seqStart as Start, c.seqLength as Length, e.exposureTime as Exposure, c.sequence as Sequence, 
                        c.modification as Modification, b.fragmentType as Fragment, a.mhp as Center, a.intensitySum as Intensity, a.retentionTime as RT, b.mhp as MHP 
                        FROM Aggregate a JOIN Component b on a.componentID=b.ID JOIN Peptide c on b.peptideID=c.ID 
                        JOIN Protein d on c.proteinID=d.ID JOIN RawFiles e on a.rawFileID=e.ID JOIN States f on e.stateID=f.ID
                        WHERE CENTER IS NOT NULL AND CENTER != 0 AND PepID = '''+str(peptide['PepID']))
                    headers = [n[0] for n in c.description]
                    pep_list = [dict(zip(headers, n)) for n in c.fetchall()]
                for state in list(set([i['State'] for i in pep_list])):
                    pep_state_list = [i for i in pep_list if i['State']==state]
                    for exposure in list(set([i['Exposure'] for i in pep_state_list])):
                        replicate_list = [i for i in pep_state_list if i['Exposure']==exposure]
                        if len(replicate_list) > 0:
                            mhp = replicate_list[0]['MHP']
                            center_list = [float(i['Center']) for i in replicate_list]
                            intensity_list = [float(i['Intensity']) for i in replicate_list]
                            isum = sum(intensity_list)
                            weighted_stats = DescrStatsW(center_list, weights=intensity_list, ddof=0)
                            center = float(weighted_stats.mean)
                            centersd = float(weighted_stats.std)
                            rt_list = [float(i['RT']) for i in replicate_list]
                            weighted_stats = DescrStatsW(rt_list, weights=intensity_list, ddof=0)
                            rt = float(weighted_stats.mean)
                            rtsd = float(weighted_stats.std)
                            sequence = replicate_list[0]['Sequence']
                            up_count = 0
                            for char in list(sequence)[1:]:
                                if char != "P":
                                    up_count += 1
                            maxuptake = int(up_count)
                            master_csv.append({'PepID': replicate_list[0]['PepID'],
                                               'Protein': replicate_list[0]['Protein'],
                                               'Start': replicate_list[0]['Start'] + 1,
                                               'End': replicate_list[0]['Start'] + replicate_list[0]['Length'],
                                               'size': replicate_list[0]['Length'],
                                               'Sequence': replicate_list[0]['Sequence'],
                                               'Modification': replicate_list[0]['Modification'],
                                               'Fragment': '',
                                               'MaxUptake': maxuptake,
                                               'MHP': mhp,
                                               'State': replicate_list[0]['State'],
                                               'Exposure': float(replicate_list[0]['Exposure']),
                                               'Center': center,
                                               'Center SD': centersd,
                                               'Uptake': 0,
                                               'Uptake SD': 0,
                                               'RT': rt,
                                               'RT SD': rtsd})
                            if replicate_list[0]['Exposure'] == 0:
                                if replicate_list[0]['State'] in zeroes.keys():
                                    if replicate_list[0]['PepID'] not in zeroes[replicate_list[0]['State']].keys():
                                        zeroes[replicate_list[0]['State']][replicate_list[0]['PepID']] = \
                                            {'State': replicate_list[0]['State'],
                                             'Center': center,
                                             'Center SD': centersd}
                                else:
                                    zeroes[replicate_list[0]['State']] = {}
                                    zeroes[replicate_list[0]['State']][replicate_list[0]['PepID']] = \
                                        {'State': replicate_list[0]['State'],
                                         'Center': center,
                                         'Center SD': centersd}
                n = n + 1
            for row in master_csv:
                if row['State'] in zeroes.keys():
                    if row['PepID'] in zeroes[row['State']].keys():
                        row['Uptake'] = row['Center'] - zeroes[row['State']][row['PepID']]['Center']
                        row['Uptake SD'] = sqrt(
                            row['Center SD'] ** 2 + zeroes[row['State']][row['PepID']]['Center SD'] ** 2)
                    else:
                        row['Uptake'] = 0
                        row['Uptake SD'] = 0
            self.master_csv = master_csv

            #
            # for state in states:
            #     for exposure in exposures:
            #         for peptide in self.peptides:
            #             peptide_list = [row for row in self.aggregate if
            #                             row['PepID'] == peptide['PepID'] and row['State'] == state and float(row[
            #                                 'Exposure']) == exposure and row['Center'] is not None]
            #             if len(peptide_list) > 0:
            #                 mhp = peptide_list[0]['MHP']
            #                 center_list = [float(i['Center']) for i in peptide_list]
            #                 intensity_list = [float(i['Intensity']) for i in peptide_list]
            #                 isum = sum(intensity_list)
            #                 weighted_stats = DescrStatsW(center_list, weights=intensity_list, ddof=0)
            #                 center = float(weighted_stats.mean)
            #                 centersd = float(weighted_stats.std)
            #                 rt_list = [float(i['RT']) for i in peptide_list]
            #                 weighted_stats = DescrStatsW(rt_list, weights=intensity_list, ddof=0)
            #                 rt = float(weighted_stats.mean)
            #                 rtsd = float(weighted_stats.std)
            #                 sequence = peptide_list[0]['Sequence']
            #                 up_count = 0
            #                 for char in list(sequence)[1:]:
            #                     if char != "P":
            #                         up_count += 1
            #                 maxuptake = int(up_count)
            #                 master_csv.append({'PepID': peptide_list[0]['PepID'],
            #                                    'Protein': peptide_list[0]['Protein'],
            #                                    'Start': peptide_list[0]['Start'] + 1,
            #                                    'End': peptide_list[0]['Start'] + peptide_list[0]['Length'],
            #                                    'size': peptide_list[0]['Length'],
            #                                    'Sequence': peptide_list[0]['Sequence'],
            #                                    'Modification': peptide_list[0]['Modification'],
            #                                    'Fragment': '',
            #                                    'MaxUptake': maxuptake,
            #                                    'MHP': mhp,
            #                                    'State': peptide_list[0]['State'],
            #                                    'Exposure': float(peptide_list[0]['Exposure']),
            #                                    'Center': center,
            #                                    'Center SD': centersd,
            #                                    'Uptake': 0,
            #                                    'Uptake SD': 0,
            #                                    'RT': rt,
            #                                    'RT SD': rtsd})
            #                 if peptide_list[0]['Exposure'] == 0:
            #                     if peptide_list[0]['State'] in zeroes.keys():
            #                         if peptide_list[0]['PepID'] not in zeroes[peptide_list[0]['State']].keys():
            #                             zeroes[peptide_list[0]['State']][peptide_list[0]['PepID']] = \
            #                                 {'State': peptide_list[0]['State'],
            #                                  'Center': center,
            #                                  'Center SD': centersd}
            #                     else:
            #                         zeroes[peptide_list[0]['State']] = {}
            #                         zeroes[peptide_list[0]['State']][peptide_list[0]['PepID']] = \
            #                             {'State': peptide_list[0]['State'],
            #                              'Center': center,
            #                              'Center SD': centersd}
            #                 n = n + 1
            #
            # for row in master_csv:
            #     if row['State'] in zeroes.keys():
            #         if row['PepID'] in zeroes[row['State']].keys():
            #             row['Uptake'] = row['Center'] - zeroes[row['State']][row['PepID']]['Center']
            #             row['Uptake SD'] = sqrt(
            #                 row['Center SD'] ** 2 + zeroes[row['State']][row['PepID']]['Center SD'] ** 2)
            #         else:
            #             row['Uptake'] = 0
            #             row['Uptake SD'] = 0
            # self.master_csv = master_csv

        else:
            exp_factors_set = set()
            self.pep_corrected = False
            with open(infile, 'r') as csvfile:
                reader = csv.DictReader(csvfile)
                timepoints = {}
                for row in reader:
                    if timepoints == {}:
                        for i in row.keys():
                            if i[0:7] == 'Uptake_':
                                if i[-1] == 'h':
                                    timepoints[i] = (float(i[7:-1]) * 60)
                                elif i[-1] == 'm':
                                    timepoints[i] = (float(i[7:-1]))
                                elif i[-1] == 's':
                                    timepoints[i] = (float(i[7:-1]) / 60)
                                else:
                                    print('Unrecognized Character')
                    for item in timepoints.keys():
                        row['Exposure'] = timepoints[item]
                        row['Uptake'] = row[item]
                        row['Uptake SD'] = row[str('UptakeSD_' + item[7:])]
                        self.master_csv.append(row)

    # Get protein sequence and organize master_csv into data_nest
    def organizeData(self):
        proteins_set = set()
        states_set = set()
        exposures_set = set()
        for row in self.master_csv:
            proteins_set.add(str(row['Protein']))
            exposures_set.add(float(row['Exposure']))
            states_set.add(row['State'])
            row['Exposure'] = float(row['Exposure'])
        self.proteins = sorted(list(proteins_set))
        self.states = sorted(list(states_set))
        self.exposures = sorted(list(exposures_set))
        # Make Protein-separated peptide list

        protein_array = {}
        for protein in self.proteins:
            protein_array[protein] = []
        for row in self.master_csv:
            if row["Protein"] in self.proteins:
                protein_array[row["Protein"]].append(row)

        # Make residue list
        sequences_list = {}
        for protein in self.proteins:
            self.residues[protein] = []
            sequences_list[protein] = []
            for row in protein_array[protein]:
                for n in range(row['Start'], row['End'] + 1):
                    if n not in self.residues[protein]:
                        self.residues[protein].append(n)
                self.residues[protein].sort()
            if min(self.residues[protein])<1:
                sequences_list[protein] = ['X'] * (abs(min(self.residues[protein]))+max(self.residues[protein])+1)
                min_val = (abs(min(self.residues[protein])))
            else:
                sequences_list[protein] = ['X'] * (max(self.residues[protein]))
                min_val = -1
            for row in protein_array[protein]:
                for n in range(row['Start'], row['End'] + 1):
                    aa_old = sequences_list[protein][n + min_val]
                    aa_new = row['Sequence'][n - row['Start']]
                    if aa_old != aa_new:
                        if aa_old == 'X':
                            sequences_list[protein][n + min_val] = aa_new
                        else:
                            print('Sequence Mismatch: AA:', n, ' Old:', aa_old, ' New:', aa_new)
            self.sequences[protein] = ''.join(sequences_list[protein])

        # Separate proteins in sorted_protein_array by time and state
        # (Dictionary of Dictionaries of Dictionaries of Lists of Dictionaries)

        for x in self.proteins:
            # Sort protein_array lists by size of peptide
            sorted_protein_array = sorted(protein_array[x], key=itemgetter('size'))
            self.data_nest[x] = {}
            for y in self.states:
                self.data_nest[x][y] = {}
                # List of peptides in a state for a specific protein
                state_array = [i for i in sorted_protein_array if i['State'] == y]
                for z in self.exposures:
                    # Dictionary of exposures with values as lists or corresponding peprides
                    self.data_nest[x][y][z] = [i for i in state_array if (float(i['Exposure']) == z) and (i['Uptake'] != nan)]

    # Delete empty states ***PROBLEMATIC: 0 may not always be performed
    def cleanFiles(self):
        '''
        Delete empty states
        :return:
        '''
        for protein in list(self.data_nest.keys()):
            for state in list(self.data_nest[protein].keys()):
                if self.data_nest[protein][state][0] == []:
                    for exposure in list(self.data_nest[protein][state].keys()):
                        for row in self.data_nest[protein][state][exposure]:
                            del self.master_csv[self.master_csv.index(row)]
                            del row
                    del self.data_nest[protein][state]

    # Sets peptide IDs for non-DynamX filetypes
    def assignIDs_CSV(self):
        self.peplist = {}
        separatedlist = []
        self.pepid=-1
        for row in sorted(sorted(sorted(self.master_csv,key=itemgetter('End')),key=itemgetter('Start')),key=itemgetter('Protein')):
            id = {'Protein':row['Protein'], 'Sequence':row['Sequence'], 'Start':row['Start'], 'End':row['End'], 'Modification':row['Modification'], 'Fragment':row['Fragment']}
            if id not in separatedlist:
                self.pepid = self.pepid + 1
                row['PepID'] = self.pepid
                self.peplist[self.pepid] = id
                separatedlist.append(id)
            else:
                row['PepID'] = separatedlist.index(id)

    # Sets peptide IDs for DynamX filetypes
    def assignIDs_DNX(self):
        self.peplist = {}
        separatedlist = []
        indices = []
        self.pepid = 0
        for row in sorted(sorted(sorted(self.master_csv, key=itemgetter('End')), key=itemgetter('Start')),
                          key=itemgetter('Protein')):
            id = {'Protein': row['Protein'], 'Sequence': row['Sequence'], 'Start': row['Start'], 'End': row['End'],
                  'Modification': row['Modification'], 'Fragment': row['Fragment']}
            if id not in separatedlist:
                if row['PepID'] is None:
                    self.pepid -= 1
                    row['PepID'] = self.pepid
                    self.peplist[self.pepid] = id
                    separatedlist.append(id)
                    indices.append(self.pepid)
                else:
                    self.peplist[row['PepID']] = id
                    separatedlist.append(id)
                    indices.append(row['PepID'])
            else:
                row['PepID'] = indices[separatedlist.index(id)]

    # Merge function, rename state
    def renameState(self,oldname,newname):
        '''
        Rename states
        :param oldname:
        :param newname:
        :return:
        '''
        for protein in self.data_nest.keys():
            self.data_nest[protein][newname] = self.data_nest[protein][oldname]
            del self.data_nest[protein][oldname]
            for exposure in self.data_nest[protein][newname].keys():
                for row in self.data_nest[protein][newname][exposure]:
                    row['State'] = newname

    # Get Ions if DynamX filetype, using specified commands
    def getDnxIons(self, commands=[], null=False):
        '''
        Returns ions based on filtering parameters
        :param commands: SQL type commands
        :type commands: list of strings
        :return: list of ions
        '''
        if str(splitext(split(self.filename)[1])[1]).lower() == '.dnx':
            if null:
                command = '''SELECT mz, 
                                    retentionTime as RT, 
                                    drift as Drift, 
                                    rawFileID as RawID, 
                                    intensity as Int 
                            FROM ions 
                            WHERE clusterID IS NULL'''
                if len(commands) > 0:
                    command = command + ' AND '
                    n = 1
                    for item in commands:
                        command = command + item
                        if len(commands) > 1 and n < len(commands):
                            command = command + ' AND '
                        n = n + 1
                with sqlite3.connect(self.filename) as conn:
                    conn.row_factory = sqlite3.Row
                    c = conn.cursor()
                    c.execute(command)
                    headers = [n[0] for n in c.description]
                    ions = [dict(zip(headers, n)) for n in c.fetchall()]
                conn.close()
                return ions
            else:
                command = '''SELECT a.mz as mz, 
                                    a.mzSearchError as searchPPM, 
                                    a.intensity as Int, 
                                    a.retentionTime as RT, 
                                    a.drift as Drift, 
                                    g.ID as RawID, 
                                    h.name as State, 
                                    f.name as Protein, 
                                    g.exposureTime as Exposure, 
                                    b.chargeState as z, 
                                    e.sequence as Sequence, 
                                    e.seqStart as Start, 
                                    e.seqLength as Length, 
                                    d.peptideID as PepID 
                            FROM ions a 
                            JOIN cluster b on a.clusterID = b.ID 
                            JOIN aggregate c on b.aggregateID = c.ID 
                            JOIN rawfiles r on c.rawFileID = r.ID 
                            JOIN component d on c.componentID = d.ID 
                            JOIN peptide e on d.peptideID = e.ID 
                            JOIN protein f on e.proteinID = f.ID 
                            JOIN rawfiles g on a.rawFileID=g.ID 
                            JOIN states h on g.stateID=h.ID'''
                if len(commands) > 0:
                    command = command + ' WHERE '
                    n = 1
                    for item in commands:
                        command = command + item
                        if len(commands) > 1 and n < len(commands):
                            command = command + ' AND '
                        n = n + 1
                with sqlite3.connect(self.filename) as conn:
                    conn.row_factory = sqlite3.Row
                    c = conn.cursor()
                    c.execute(command)
                    headers = [n[0] for n in c.description]
                    ions = [dict(zip(headers, n)) for n in c.fetchall()]
                conn.close()
                return ions
        else:
            print('File not DnX Format')
            return []

    def getDnxCluster(self,commands=[]):
        if str(splitext(split(self.filename)[1])[1]).lower() == '.dnx':
            command = '''SELECT b.rawFileID as RawID, 
                                d.ID as PepID, 
                                a.ID as ClusterID, 
                                b.ID as AggregateID, 
                                c.ID as ComponentID, 
                                e.name as Protein, 
                                c.seqStart as Start, 
                                c.seqLength as Length, 
                                c.sequence as Sequence, 
                                c.modification as Modification, 
                                c.fragmentType as Fragment,
                                c.maxD as MaxUptake, 
                                c.mhp as MHP, 
                                g.name as State, 
                                f.exposureTime as Exposure, 
                                a.chargeState as z, 
                                b.mhp as Center, 
                                b.mhpVariance as CenterSD, 
                                a.mzAverage as mzCenter,
                                a.massStart as mzCenterStart, 
                                a.massEnd as mzCenterEnd, 
                                b.retentionTime as RT, 
                                b.retentionTimeFWHM as RTFWHM, 
                                a.timeStart as RTStart, 
                                a.timeEnd as RTEnd,
                                a.driftAverage as Drift, 
                                a.driftFWHM as DriftFWHM, 
                                a.driftStart as DriftStart, 
                                a.driftEnd as DriftEnd
                    FROM Cluster a 
                    JOIN Aggregate b on a.aggregateID=b.ID 
                    JOIN Component c on b.componentID=c.ID 
                    JOIN Peptide d on c.peptideID=d.ID 
                    JOIN Protein e on d.proteinID=e.ID 
                    JOIN Rawfiles f on b.rawFileID=f.ID 
                    JOIN States g on f.stateID=g.ID 
                    WHERE a.mzAverage IS NOT NULL AND a.mzAverage != 0 AND a.scanStart != -1 AND a.scanEnd !=-1'''
            if len(commands) > 0:
                command = command + ' AND '
                n = 1
                for item in commands:
                    command = command + item
                    if len(commands) > 1 and n < len(commands):
                        command = command + ' AND '
                    n = n + 1
            with sqlite3.connect(self.filename) as conn:
                conn.row_factory = sqlite3.Row
                c = conn.cursor()
                c.execute(command)
                headers = [n[0] for n in c.description]
                cluster = [dict(zip(headers, n)) for n in c.fetchall()]
            return cluster
        else:
            print('File not DnX Format')
            return []

    # Import Cluster, Aggregate, Ions, Peptides and Rawfiles tables
    def aquireIons(self):
        with sqlite3.connect(self.filename) as conn:
            conn.row_factory = sqlite3.Row
            c = conn.cursor()
            if self.que is not None:
                self.que.put(0)

            c.execute(
                "SELECT a.ID as PepID, c.name as Protein, a.sequence as Sequence,a.modification as Modification, b.minZ as minZ, b.maxZ as maxZ FROM Peptide a JOIN Component b on a.ID = b.peptideID JOIN Protein c on a.proteinID = c.ID WHERE scoreDB IS NOT NULL")
            headers = [n[0] for n in c.description]
            peptides = [dict(zip(headers, n)) for n in c.fetchall()]

            c.execute(
                "SELECT a.ID as RawID, b.name as State, a.exposureTime as Exposure FROM Rawfiles a JOIN States b on a.stateID=b.ID")
            headers = [n[0] for n in c.description]
            rawfiles = [dict(zip(headers, n)) for n in c.fetchall()]

            if self.que is not None:
                self.que.put(25)
        conn.close()
        self.peptides = peptides
        self.rawfiles = rawfiles

    # Calculate statistical significance for entire dataset, or just specified pepid if not None
    def getSignificance(self,pepid = None):
        if self.filetype == 'DNX':
            file_nest = {}
            for exposure in self.exposures:
                file_nest[exposure] = {i: [] for i in self.states}
            for rawfile in self.rawfiles:
                file_nest[rawfile['Exposure']][rawfile['State']].append(rawfile['RawID'])
            print('Files: ' + str(file_nest))
            stats_nest = {}
            if pepid is not None:
                peplist = [i for i in self.peptides if i['PepID']==pepid]
            else:
                peplist = self.peptides
            for peptide in peplist:
                print(peptide['Sequence'])
                stats_nest[peptide['PepID']] = {i: {} for i in self.exposures}
                print(peptide['PepID'])
                with sqlite3.connect(self.filename) as conn:
                    conn.row_factory = sqlite3.Row
                    c = conn.cursor()
                    c.execute('''SELECT a.mzAverage as mzAverage, a.chargeState as z, b.rawFileID as RawID FROM Cluster a JOIN Aggregate b on a.aggregateID=b.ID JOIN Component c on 
                       b.componentID=c.ID WHERE mzAverage IS NOT NULL AND c.peptideID = ?''',(peptide['PepID'],))
                    headers = [n[0] for n in c.description]
                    pepcluster = [dict(zip(headers, n)) for n in c.fetchall()]
                conn.close()
                for exposure in file_nest.keys():
                    print(exposure)
                    stats_nest[peptide['PepID']][exposure] = {}
                    stats_nest[peptide['PepID']][exposure]['Stats'] = {}
                    for state in set(file_nest[exposure].keys()) & set(self.data_nest[peptide['Protein']].keys()):
                        stats_nest[peptide['PepID']][exposure][state] = []
                        print(state)
                        for rawfile in file_nest[exposure][state]:
                            cluster = [i for i in pepcluster if i['RawID']==rawfile]
                            for row in cluster:
                                if (row['mzAverage'] is not None) and (row['mzAverage'] != 0):
                                    mz0 = [i['Center'] for i in self.data_nest[peptide['Protein']][state][0] if
                                           i['PepID'] == peptide['PepID']]
                                    if len(mz0) != 0:
                                        uptake = (float(row['mzAverage']) * float(row['z']) - (
                                                float(row['z']) - 1) * 1.007825) - mz0[0]
                                        stats_nest[peptide['PepID']][exposure][state].append(uptake)
                                    else:
                                        print('Error')

                        if stats_nest[peptide['PepID']][exposure][state] == []:
                            del stats_nest[peptide['PepID']][exposure][state]

                    full_states = [i for i in stats_nest[peptide['PepID']][exposure].keys() if i != 'Stats']
                    data_nest = [stats_nest[peptide['PepID']][exposure][state] for state in full_states]
                    if len(full_states) > 2:
                        print(full_states)
                        F_levene, p_levene = levene(*data_nest, center='mean')
                        F_anova, p_anova = f_oneway(*data_nest)
                        data_list = []
                        for state in full_states:
                            data_list.extend([(state, i) for i in stats_nest[peptide['PepID']][exposure][state]])
                        groups = array([i[0] for i in data_list])
                        data = array([i[1] for i in data_list])
                        mc = MultiComparison(data, groups)
                        result = mc.tukeyhsd()
                        headers = result._results_table.data[0]
                        result_dict = [dict(zip(headers, n)) for n in result._results_table.data[1:]]
                        for r in result_dict:
                            r['p-value'] = r['p-adj']
                        stats_nest[peptide['PepID']][exposure]['Stats']['Levene'] = {'F': F_levene, 'p': p_levene}
                        stats_nest[peptide['PepID']][exposure]['Stats']['Anova'] = {'F': F_anova, 'p': p_anova}
                        stats_nest[peptide['PepID']][exposure]['Stats']['Tukey'] = result_dict
                        stats_nest[peptide['PepID']][exposure]['Stats']['Tukey_obj'] = result
                    elif len(full_states) == 2:
                        ttest = ttest_ind(data_nest[0],data_nest[1])
                        means = mean(data_nest[0])-mean(data_nest[1])
                        reject = ttest[1]<0.05

                        k = len(full_states)
                        N = sum([len(i) for i in data_nest])
                        qcrit = qsturng(0.95, N, N - k)
                        msw = sum([(std(i) ** 2) * len(i) for i in data_nest]) / (N - k)
                        sci = simultaneous_ci(qcrit, msw, [len(i) for i in data_nest])
                        stats_nest[peptide['PepID']][exposure]['Stats']['CI'] = {i: sci[index] for index, i in
                                                                                 enumerate(full_states)}
                        stats_nest[peptide['PepID']][exposure]['Stats']['Ttest'] = {'group1':full_states[0],
                                                                                     'group2':full_states[1],
                                                                                     'meandiff':means,
                                                                                     'reject':str(reject),
                                                                                     'p-value':ttest[1]}
                    else:
                        stats_nest[peptide['PepID']][exposure]['Stats'] = {}
            self.stats_nest = stats_nest
        else:
            if self.corrected == 'Yes':
                uptakekey = 'Uptake_corr'
                sdkey = 'Uptake_SD_corr'
            else:
                uptakekey = 'Uptake'
                sdkey = 'Uptake SD'
            # not finished
            stats_nest = {}
            n = 1
            if pepid is not None:
                peplist = {pepid:self.peplist[pepid]}
            else:
                peplist = self.peplist
            for (index,peptide) in [(i,peplist[i]) for i in peplist.keys()]:
                stats_nest[index] = {i: {} for i in self.exposures}
                for exposure in [e for e in self.exposures if e != 0]:
                    stats_nest[index][exposure]['Stats'] = {}
                    for state in self.data_nest[peptide['Protein']].keys():
                        uptake = [i for i in self.data_nest[peptide['Protein']][state][exposure] if i['PepID']==index]
                        if len(uptake) > 0:
                            stats_nest[index][exposure][state] = (float(uptake[0][uptakekey]), float(uptake[0][sdkey]),self.replicates[state][exposure])
                    full_states = [i for i in stats_nest[index][exposure].keys() if (i != 'Stats') and (stats_nest[index][exposure][i][2]>1)]
                    data_l = [stats_nest[index][exposure][state] for state in full_states]
                    
                    if len(full_states) > 2:
                        k = len(full_states)
                        N = sum([self.replicates[i][exposure] for i in self.replicates.keys()])
                        df_between = k - 1
                        df_within = N - k
                        grand_mean = sum([i[0] for i in data_l]) / k
                        variance_means = sum([(i[0] - grand_mean) ** 2 for i in data_l]) / df_between
                        meansquare_between = variance_means * mean([self.replicates[i][exposure] for i in self.replicates.keys()])
                        meansquare_within = sum([i[1] for i in data_l]) / k
                        F = meansquare_between / meansquare_within
                        crit_anova = f.ppf(q=0.95, dfn=df_between, dfd=df_within)
                        p = f.sf(F, dfn=df_between, dfd=df_within)
                        stats_nest[index][exposure]['Stats']['Anova'] = {'F': F, 'p': p, 'Critical': crit_anova}

                        state_combos = list(combinations(full_states, 2))
                        qcrit = qsturng(0.95, N, N - k)
                        msw = sum([(i[1] ** 2) * i[2] for i in data_l]) / (N-k)
                        mean_diffs = [stats_nest[index][exposure][i[0]][0] - stats_nest[index][exposure][i[1]][0] for
                            i in state_combos]
                        std_pairs = [sqrt(msw*(1/self.replicates[i[0]][exposure]+1/self.replicates[i[1]][exposure])/2) for i in state_combos]
                        p_values = psturng(abs(mean_diffs)/std_pairs,k,N-k)
                        sci = simultaneous_ci(qcrit, msw, [self.replicates[i][exposure] for i in full_states])
                        result = []
                        for index2,group in enumerate(state_combos):
                            row = {}
                            row['group1'] = group[0]
                            row['group2'] = group[1]
                            row['HSD'] = std_pairs[index2] * qcrit
                            row['meandiff'] = mean_diffs[index2]
                            row['lower'] = row['meandiff'] - row['HSD']
                            row['upper'] = row['meandiff'] + row['HSD']
                            row['p-value'] = p_values[index2]
                            if float(row['p-value']) < 0.05:
                                row['reject'] = True
                            else:
                                row['reject'] = False
                            result.append(row)
                        stats_nest[index][exposure]['Stats']['CI'] = {i:sci[index] for index,i in enumerate(full_states)}
                        stats_nest[index][exposure]['Stats']['Tukey'] = result
                    elif len(full_states) == 2:
                        ttest = ttest_ind_from_stats(data_l[0][0],data_l[0][1],data_l[0][2],
                                                     data_l[1][0],data_l[1][1],data_l[1][2])
                        means = data_l[0][0] - data_l[1][0]
                        reject = ttest[1] < 0.05
                        k = len(full_states)
                        N = sum([self.replicates[i][exposure] for i in self.replicates.keys()])
                        qcrit = qsturng(0.95, N, N - k)
                        msw = sum([(i[1] ** 2) * i[2] for i in data_l]) / (N - k)
                        sci = simultaneous_ci(qcrit, msw, [self.replicates[i][exposure] for i in full_states])
                        stats_nest[index][exposure]['Stats']['CI'] = {i:sci[index] for index,i in enumerate(full_states)}
                        stats_nest[index][exposure]['Stats']['Ttest'] = {'group1': full_states[0],
                                                                                    'group2': full_states[1],
                                                                                    'meandiff': means,
                                                                                    'reject': str(reject),
                                                                                    'p-value': ttest[1]}
                    else:
                        stats_nest[index][exposure]['Stats'] = {}
                n = n + 1
            self.stats_nest = stats_nest

    # Find outlier peak assignments in data from DNX project files
    def findOutliers(self, pep = None):
        if pep is not None:
            peplist = [i for i in self.peptides if i['PepID']==pep]
        else:
            peplist = self.peptides
        file_nest = {}

        for peptide in peplist:
            pepions = self.getDnxIons(['a.clusterID IS NOT NULL','PepID = '+str(peptide['PepID'])])
            charges = [z for z in range(int(peptide['minZ']), int(peptide['maxZ']) + 1)]
            file_nest[peptide['PepID']] = {'_STATS_': {}, 'States': {}, 'Exposures': {}}
            pep_states = set([str(i['State']) for i in pepions])
            for state in pep_states:
                state_exposures = [float(i['Exposure']) for i in self.rawfiles if str(i['State']) == state]
                file_nest[peptide['PepID']]['States'][state] = {'_STATS_': {}}
                for exp in state_exposures:
                    file_nest[peptide['PepID']]['States'][state][exp] = {'_STATS_': {}}
                    exposure_raws = [i for i in self.rawfiles if
                                     (i['State'] == str(state)) and (float(i['Exposure']) == exp)]
                    for raw in exposure_raws:
                        file_nest[peptide['PepID']]['States'][state][exp][raw['RawID']] = {'_STATS_': {}}
                        for z in charges:
                            file_nest[peptide['PepID']]['States'][state][exp][raw['RawID']][z] = {}
            for exp in self.exposures:
                file_nest[peptide['PepID']]['Exposures'][exp] = {'_STATS_': {}}
            id = peptide['PepID']
            states = list(file_nest[peptide['PepID']]['States'].keys())
            for state in states:
                stateions = [i for i in pepions if (str(i['State']) == str(state))]
                exps = [i for i in file_nest[peptide['PepID']]['States'][state].keys() if i != '_STATS_']
                for exp in exps:
                    expions = [i for i in stateions if (float(i['Exposure']) == float(exp))]
                    raws = [i for i in file_nest[peptide['PepID']]['States'][state][exp].keys() if i != '_STATS_']
                    for raw in raws:
                        rawions = [i for i in expions if (int(i['RawID']) == int(raw))]
                        charges = [i for i in file_nest[peptide['PepID']]['States'][state][exp][raw].keys() if
                                   i != '_STATS_']
                        for z in charges:
                            zions = [i for i in rawions if (int(i['z']) == int(z))]
                            if len(zions) == 0:
                                file_nest[id]['States'][state][exp][raw][z] = {}
                            else:
                                mz_list = [i['mz'] for i in zions]
                                ppm_list = [i['searchPPM'] for i in zions]
                                rt_list = [i['RT'] for i in zions]
                                drift_list = [i['Drift'] for i in zions]
                                mz_list, ppm_list, rt_list, drift_list = [list(i) for i in
                                                                          zip(*sorted(zip(mz_list, ppm_list, rt_list,
                                                                                          drift_list)))]

                                std_rt = std(rt_list)
                                mean_rt = mean(rt_list)
                                range_rt = (mean_rt - 2 * std_rt, mean_rt + 2 * std_rt)

                                std_drift = std(drift_list)
                                mean_drift = mean(drift_list)
                                range_drift = (mean_drift - 2 * std_drift, mean_drift + 2 * std_drift)

                                file_nest[id]['States'][state][exp][raw][z] = \
                                    {'data': mz_list,
                                     'ppm': {'data': ppm_list},
                                     'rt': {'data': rt_list,
                                            'mean':mean_rt,
                                            'std':std_rt},
                                     'drift': {'data': drift_list,
                                               'mean': mean_drift,
                                               'std': std_drift}}
                        # RawFile Stats
                        for z in [i for i in file_nest[peptide['PepID']]['States'][state][exp][raw].keys() if
                                  i != '_STATS_']:
                            if list(file_nest[id]['States'][state][exp][raw][z].keys()) == []:
                                del file_nest[id]['States'][state][exp][raw][z]
                        new_charges = [z for z in file_nest[id]['States'][state][exp][raw].keys() if z != '_STATS_']
                        if len(new_charges) == 0:
                            print('Peptide not found in rawfile')
                            del file_nest[id]['States'][state][exp][raw]
                        else:
                            mz_lists = [file_nest[id]['States'][state][exp][raw][z]['data'] for z in new_charges]
                            mz_list = [item for sublist in mz_lists for item in sublist]
                            ppm_lists = [file_nest[id]['States'][state][exp][raw][z]['ppm']['data'] for z in
                                         new_charges]
                            ppm_list = [item for sublist in ppm_lists for item in sublist]
                            rt_lists = [file_nest[id]['States'][state][exp][raw][z]['rt']['data'] for z in new_charges]
                            rt_list = [item for sublist in rt_lists for item in sublist]
                            mz_list, ppm_list, rt_list = [list(i) for i in
                                                          zip(*sorted(zip(mz_list, ppm_list, rt_list)))]
                            drift_list = {}
                            for z in new_charges:
                                if z not in list(drift_list.keys()):
                                    drift_list[z] = []
                                drift_list[z].extend(
                                    file_nest[id]['States'][state][exp][raw][z]['drift']['data'])
                            std_rt = std(rt_list)
                            mean_rt = mean(rt_list)
                            range_rt = (mean_rt - 2 * std_rt, mean_rt + 2 * std_rt)
                            iqr_drift = {}
                            std_drift = {}
                            mean_drift = {}
                            range_drift = {}
                            for z in list(drift_list.keys()):

                                std_drift[z] = std(drift_list[z])
                                mean_drift[z] = mean(drift_list[z])
                                range_drift[z] = (mean_drift[z] - 2 * std_drift[z], mean_drift[z] + 2 * std_drift[z])

                            file_nest[id]['States'][state][exp][raw]['_STATS_'] = \
                                {'data': mz_list,
                                 'ppm': {'data': ppm_list},
                                 'rt': {'data': rt_list, 'mean':mean_rt,'std':std_rt},
                                 'drift': {'data': drift_list,'mean':mean_drift,'std':std_drift}}

                    # Exposure Stats
                    for raw in [i for i in file_nest[peptide['PepID']]['States'][state][exp].keys() if i != '_STATS_']:
                        if list(file_nest[id]['States'][state][exp][raw].keys()) == ['_STATS_']:
                            del file_nest[id]['States'][state][exp][raw]
                    new_raws = [i for i in file_nest[id]['States'][state][exp].keys() if i != '_STATS_']
                    if len(new_raws) == 0:
                        print('Peptide not found in exposure')
                        del file_nest[id]['States'][state][exp]
                    else:
                        mz_lists = [file_nest[id]['States'][state][exp][raw][z]['data']
                                    for raw in new_raws
                                    for z in
                                    [i for i in file_nest[id]['States'][state][exp][raw].keys() if i != '_STATS_']]
                        mz_list = [item for sublist in mz_lists for item in sublist]
                        ppm_lists = [file_nest[id]['States'][state][exp][raw][z]['ppm']['data']
                                     for raw in new_raws
                                     for z in
                                     [i for i in file_nest[id]['States'][state][exp][raw].keys() if i != '_STATS_']]
                        ppm_list = [item for sublist in ppm_lists for item in sublist]
                        rt_lists = [file_nest[id]['States'][state][exp][raw][z]['rt']['data']
                                    for raw in new_raws
                                    for z in
                                    [i for i in file_nest[id]['States'][state][exp][raw].keys() if i != '_STATS_']]
                        rt_list = [item for sublist in rt_lists for item in sublist]
                        mz_list, ppm_list, rt_list = [list(i) for i in
                                                      zip(*sorted(zip(mz_list, ppm_list, rt_list)))]
                        drift_list = {}
                        for raw in new_raws:
                            for z in [i for i in file_nest[id]['States'][state][exp][raw].keys() if i != '_STATS_']:
                                if z not in list(drift_list.keys()):
                                    drift_list[z] = []
                                drift_list[z].extend(file_nest[id]['States'][state][exp][raw][z]['drift']['data'])

                        std_rt = std(rt_list)
                        mean_rt = mean(rt_list)
                        range_rt = (mean_rt - 2 * std_rt, mean_rt + 2 * std_rt)
                        iqr_drift = {}
                        std_drift = {}
                        mean_drift = {}
                        range_drift = {}
                        for z in list(drift_list.keys()):

                            std_drift[z] = std(drift_list[z])
                            mean_drift[z] = mean(drift_list[z])
                            range_drift[z] = (mean_drift[z] - 2 * std_drift[z], mean_drift[z] + 2 * std_drift[z])

                        file_nest[id]['States'][state][exp]['_STATS_'] = \
                            {'data': mz_list,
                             'ppm': {'data': ppm_list},
                             'rt': {'data': rt_list, 'mean': mean_rt, 'std': std_rt},
                             'drift': {'data': drift_list, 'mean': mean_drift, 'std': std_drift}}
                # States Stats
                for exp in [i for i in file_nest[peptide['PepID']]['States'][state].keys() if i != '_STATS_']:
                    if list(file_nest[id]['States'][state][exp].keys()) == ['_STATS_']:
                        del file_nest[id]['States'][state][exp]
                new_exps = [i for i in file_nest[id]['States'][state].keys() if i != '_STATS_']
                if len(new_exps) == 0:
                    print('Peptide not found in state')
                    del file_nest[id]['States'][state]
                else:
                    mz_lists = [file_nest[id]['States'][state][exp][raw][z]['data']
                                for exp in new_exps
                                for raw in [i for i in file_nest[id]['States'][state][exp].keys() if i != '_STATS_']
                                for z in [i for i in file_nest[id]['States'][state][exp][raw].keys() if i != '_STATS_']]
                    mz_list = [item for sublist in mz_lists for item in sublist]
                    ppm_lists = [file_nest[id]['States'][state][exp][raw][z]['ppm']['data']
                                 for exp in new_exps
                                 for raw in [i for i in file_nest[id]['States'][state][exp].keys() if i != '_STATS_']
                                 for z in
                                 [i for i in file_nest[id]['States'][state][exp][raw].keys() if i != '_STATS_']]
                    ppm_list = [item for sublist in ppm_lists for item in sublist]
                    rt_lists = [file_nest[id]['States'][state][exp][raw][z]['rt']['data']
                                for exp in new_exps
                                for raw in [i for i in file_nest[id]['States'][state][exp].keys() if i != '_STATS_']
                                for z in [i for i in file_nest[id]['States'][state][exp][raw].keys() if i != '_STATS_']]
                    rt_list = [item for sublist in rt_lists for item in sublist]
                    mz_list, ppm_list, rt_list = [list(i) for i in
                                                  zip(*sorted(zip(mz_list, ppm_list, rt_list)))]
                    drift_list = {}

                    for exp in new_exps:
                        for raw in [i for i in file_nest[id]['States'][state][exp].keys() if i != '_STATS_']:
                            for z in [i for i in file_nest[id]['States'][state][exp][raw].keys() if i != '_STATS_']:
                                if z not in list(drift_list.keys()):
                                    drift_list[z] = []
                                drift_list[z].extend(file_nest[id]['States'][state][exp][raw][z]['drift']['data'])

                    std_rt = std(rt_list)
                    mean_rt = mean(rt_list)
                    range_rt = (mean_rt - 2 * std_rt, mean_rt + 2 * std_rt)
                    iqr_drift = {}
                    std_drift = {}
                    mean_drift = {}
                    range_drift = {}
                    for z in list(drift_list.keys()):

                        std_drift[z] = std(drift_list[z])
                        mean_drift[z] = mean(drift_list[z])
                        range_drift[z] = (mean_drift[z] - 2 * std_drift[z], mean_drift[z] + 2 * std_drift[z])

                    file_nest[id]['States'][state]['_STATS_'] = \
                        {'data': mz_list,
                         'ppm': {'data': ppm_list},
                         'rt': {'data': rt_list, 'mean': mean_rt, 'std': std_rt},
                         'drift': {'data': drift_list, 'mean': mean_drift, 'std': std_drift}}
            # Peptide Stats
            for state in file_nest[peptide['PepID']]['States'].keys():
                if list(file_nest[id]['States'][state].keys()) == []:
                    del file_nest[id]['States'][state]
            new_states = list(file_nest[id]['States'].keys())
            if len(new_states) == 0:
                print('Peptide not found')
                del file_nest[id]
            else:
                mz_lists = [file_nest[id]['States'][state][exp][raw][z]['data']
                            for state in new_states
                            for exp in [i for i in file_nest[id]['States'][state].keys() if i != '_STATS_']
                            for raw in [i for i in file_nest[id]['States'][state][exp].keys() if i != '_STATS_']
                            for z in [i for i in file_nest[id]['States'][state][exp][raw].keys() if i != '_STATS_']]
                mz_list = [item for sublist in mz_lists for item in sublist]
                ppm_lists = [file_nest[id]['States'][state][exp][raw][z]['ppm']['data']
                             for state in new_states
                             for exp in [i for i in file_nest[id]['States'][state].keys() if i != '_STATS_']
                             for raw in [i for i in file_nest[id]['States'][state][exp].keys() if i != '_STATS_']
                             for z in [i for i in file_nest[id]['States'][state][exp][raw].keys() if i != '_STATS_']]
                ppm_list = [item for sublist in ppm_lists for item in sublist]
                rt_lists = [file_nest[id]['States'][state][exp][raw][z]['rt']['data']
                            for state in new_states
                            for exp in [i for i in file_nest[id]['States'][state].keys() if i != '_STATS_']
                            for raw in [i for i in file_nest[id]['States'][state][exp].keys() if i != '_STATS_']
                            for z in [i for i in file_nest[id]['States'][state][exp][raw].keys() if i != '_STATS_']]
                rt_list = [item for sublist in rt_lists for item in sublist]
                mz_list, ppm_list, rt_list = [list(i) for i in
                                              zip(*sorted(zip(mz_list, ppm_list, rt_list)))]
                drift_list = {}
                for state in new_states:
                    for exp in [i for i in file_nest[id]['States'][state].keys() if i != '_STATS_']:
                        for raw in [i for i in file_nest[id]['States'][state][exp].keys() if i != '_STATS_']:
                            for z in [i for i in file_nest[id]['States'][state][exp][raw].keys() if i != '_STATS_']:
                                if z not in list(drift_list.keys()):
                                    drift_list[z] = []
                                drift_list[z].extend(file_nest[id]['States'][state][exp][raw][z]['drift']['data'])

                std_rt = std(rt_list)
                mean_rt = mean(rt_list)
                range_rt = (mean_rt - 2 * std_rt, mean_rt + 2 * std_rt)
                iqr_drift = {}
                std_drift = {}
                mean_drift = {}
                range_drift = {}
                for z in list(drift_list.keys()):

                    std_drift[z] = std(drift_list[z])
                    mean_drift[z] = mean(drift_list[z])
                    range_drift[z] = (mean_drift[z] - 2 * std_drift[z], mean_drift[z] + 2 * std_drift[z])

                file_nest[id]['_STATS_'] = \
                    {'data': mz_list,
                     'ppm': {'data': ppm_list},
                     'rt': {'data': rt_list, 'mean': mean_rt, 'std': std_rt},
                     'drift': {'data': drift_list, 'mean': mean_drift, 'std': std_drift}}
        self.file_nest = file_nest
        return file_nest
        print('Acquired')

    # Calculate summary statistics for assigned ions
    def ionstatsum(self):
        sumrt = 0
        totalrt = 0
        sumdrift = {i:0 for i in range(1, 11)}
        totaldrift = {i:0 for i in range(1, 11)}
        for id in self.file_nest:
            mean = self.file_nest[id]['_STATS_']['rt']['mean']
            sumrt += len([i for i in self.file_nest[id]['_STATS_']['rt']['data'] if abs(float(i) - mean) / mean < 0.05])
            totalrt += len(self.file_nest[id]['_STATS_']['rt']['data'])
            for z in self.file_nest[id]['_STATS_']['drift']['data'].keys():
                meanz = self.file_nest[id]['_STATS_']['drift']['mean'][z]
                sumdrift[z] += len(
                    [i for i in self.file_nest[id]['_STATS_']['drift']['data'][z] if abs(float(i) - meanz) / meanz < 0.05])
                totaldrift[z] += len(self.file_nest[id]['_STATS_']['drift']['data'][z])

    # Merge Files according to given instructions
    def mergeFiles(self, dictlist, mergeeach):
        self.data_nest = {}
        self.master_csv = []
        for state in dictlist.keys():
            for protein in dictlist[state].keys():
                for exposure in dictlist[state][protein].keys():
                    for item in dictlist[state][protein][exposure]:
                        self.master_csv.extend(mergeeach[item].data_nest[protein][state][exposure])
        exp_factors_set = set()
        for row in self.master_csv:
            if 'ExposureCorrection' in row.keys():
                exp_factors_set.add((row['Exposure'], row['ExposureCorrection']))
            if 'PeptideCorrection' in row.keys():
                if float(row['PeptideCorrection']) != 1:
                    self.pep_corrected = True
                    self.corrected = 'Yes'
        if any(item != 1 for item in [i[1] for i in exp_factors_set]):
            self.exp_corrected = True
            self.corrected = 'Yes'
            print('any yes!')
        exp_factors_list = [item[0] for item in list(exp_factors_set)]
        if len(set([x for x in exp_factors_list if exp_factors_list.count(x) > 1])) > 0:
            self.merge_and_corrected = True
        for i in sorted(list(exp_factors_set)):
            self.cor_factor[float(i[0])] = float(i[1])
        self.organizeData()
        self.cleanFiles()

    # Populate sequence of residues with the best peptide information at each location
    def assignData(self, protein, state1, exposure, state2=None):
        data_nest = self.data_nest
        residues = self.residues
        seq_nest = []
        for i in residues[protein]:
            seq_dict = {'pos': i, 'pepmin': 0, 'pepmax': 0, 'size': 0, 'assigned': 0, 'fragment':''}
            seq_nest.append(seq_dict)
        peptide_list = data_nest[protein][state1][exposure]
        if state2 is not None:
            merge_list = []
            for i in data_nest[protein][state1][exposure]:
                for j in data_nest[protein][state2][exposure]:
                    if i["Start"] == j["Start"] and i["End"] == j["End"]:
                        merge_list.append(i)
            peptide_list = merge_list
        for n in seq_nest:
            peptidelist = [i for i in peptide_list if (i['Start'] < n['pos'] <= i['End'] and i['Fragment'] != 'C')]
            peptidelist.extend([i for i in peptide_list if (i['Start'] <= n['pos'] <= i['End'] and i['Fragment'] == 'C')])
            for peptide in peptidelist:
                if peptide['size'] < n['size'] or n['assigned'] == 0:
                    n['pepmin'] = peptide['Start']
                    n['pepmax'] = peptide['End']
                    n['size'] = peptide['size']
                    n['assigned'] = 1
                    n['fragment'] = peptide['Fragment']

        # Count Assigned Residues
        count = 0
        for i in seq_nest:
            if i['assigned'] == 1:
                count += 1
        self.seq_nest = seq_nest

    # Back Exchange Correction
    def backCorrect(self, cor_dict, global_bx=None, timepoint=None):
        self.cor_factor = cor_dict
        if self.recombined:
            for item in [i for i in self.master_csv if i['Fragment'] != '']:
                self.master_csv.remove(item)
            self.organizeData()
            self.cleanFiles()
        for protein in self.proteins:
            for state in self.data_nest[protein].keys():
                for exposure in self.data_nest[protein][state].keys():
                    for row in self.data_nest[protein][state][exposure]:
                        if row['Fragment'] != '':
                            print('Found Fragment!'+str(protein)+" "+str(state)+" "+str(exposure))
                        else:
                            if timepoint:
                                self.pep_corrected = True
                                timepoint_factor = [(float(i['Uptake'])/float(i['MaxUptake'])) for i in
                                                    self.data_nest[protein][state][float(timepoint)] if
                                                    (row['Start'] == i['Start']) and (row['End'] == i['End']) and
                                                    (row['Sequence'] == i['Sequence']) and
                                                    (row['Fragment'] == i['Fragment'])]
                                if timepoint_factor != []:
                                    timepoint_factor = float(timepoint_factor[0])
                                    if timepoint_factor == 0:
                                        timepoint_factor = 1
                                    factor = (float(cor_dict[exposure])) * timepoint_factor
                                    row['PeptideCorrection'] = timepoint_factor
                            elif global_bx is not None:
                                factor = (float(cor_dict[exposure]))*(1-(float(global_bx)/100))
                                row['GlobalCorrection'] = 1-(float(global_bx)/100)
                            else:
                                factor = float(cor_dict[exposure])
                                row['PeptideCorrection'] = 1
                            row['Uptake_corr'] = float(row['Uptake']) / factor
                            row['Uptake_SD_corr'] = float(row['Uptake SD']) / factor
                            row['ExposureCorrection'] = float(cor_dict[exposure])
                            if float(row['MaxUptake']) == 0:
                                row['Percent_Uptake_corr'] = 0
                            else:
                                row['Percent_Uptake_corr'] = float(row['Uptake_corr']) / float(row['MaxUptake']) * 100

        self.corrected = 'Yes'
        if any(cor_dict[i] != 1 for i in cor_dict.keys()):
            self.exp_corrected = True
        if self.recombined:
            self.recombined = False
            self.segmentPeptides()
            self.applySegmentation()

    # Overlapping Peptide Segmentation (OPS) initial calculations
    def segmentPeptides(self):
        res = 6
        recombine_dict = {}
        for protein in self.proteins:
            recombine_dict[protein] = {}
            for state in self.data_nest[protein].keys():
                recombine_dict[protein][state] = {}
                for exposure in self.data_nest[protein][state].keys():
                    if self.corrected == 'Yes':
                        datalist = self.data_nest[protein][state][exposure]
                        if len(datalist)>0:
                            expf = datalist[0]['ExposureCorrection']
                            pepf = 1
                        else:
                            print('Empty Exposure:'+str(protein)+' '+str(state)+' '+str(exposure))
                    new_dict = {}
                    recombine_dict[protein][state][exposure] = []
                    data = self.data_nest[protein][state][exposure]

                    # Find matching N-termini
                    # Make dictionary with C values for N keys
                    x_dict = {}
                    for x in data:
                        if x['Start'] not in x_dict.keys():
                            x_dict[x['Start']] = {x['End']: x}
                        else:
                            x_dict[x['Start']][x['End']] = x
                    # Iterate C values in N keys to create new peptides
                    for x in x_dict.keys():
                        srt_y = sorted(x_dict[x])
                        if len(srt_y) > 1:
                            for i in range(0, len(srt_y) - 1): # Changed len(srt_y) - 1 to len(srt_y) to include last item. 20180909
                                for j in srt_y[i + 1:]:
                                    uptake = round(float(x_dict[x][j]['Uptake']), res) - round(float(x_dict[x][srt_y[i]]['Uptake']), res)
                                    sd = sqrt(round(float(x_dict[x][j]['Uptake SD']), res) ** 2 + round(float(x_dict[x][srt_y[i]]['Uptake SD']), res) ** 2)
                                    seq = ''.join(x_dict[x][j]['Sequence'][-(j - srt_y[i]):])
                                    if self.corrected == 'Yes':
                                        if 'Uptake_corr' not in x_dict[x][j] or 'Uptake_corr' not in x_dict[x][srt_y[i]]:
                                            print('Error!')
                                            continue
                                        uptake_corr = round(float(x_dict[x][j]['Uptake_corr']), res) - round(float(x_dict[x][srt_y[i]]['Uptake_corr']), res)
                                        sd_corr = sqrt(round(float(x_dict[x][j]['Uptake_SD_corr']), res)**2 + round(float(x_dict[x][srt_y[i]]['Uptake_SD_corr']), res)**2)

                                        if (srt_y[i] + 1) not in new_dict.keys():
                                            new_dict[(srt_y[i] + 1)] = {j: [{'Uptake': uptake,
                                                                             'Uptake_corr': uptake_corr,
                                                                             'Uptake SD': sd, 'Uptake_SD_corr': sd_corr,
                                                                             'Sequence': seq, 'Fragment':'C'}]}
                                        else:
                                            if int(j) not in new_dict[srt_y[i] + 1].keys():
                                                new_dict[srt_y[i] + 1][j] = [{'Uptake': uptake,
                                                                              'Uptake_corr': uptake_corr,
                                                                              'Uptake SD': sd,
                                                                              'Uptake_SD_corr': sd_corr,
                                                                              'Sequence': seq,
                                                                              'Fragment':'C'}]
                                            else:
                                                new_dict[srt_y[i] + 1][j].extend([{'Uptake': uptake,
                                                                                   'Uptake_corr': uptake_corr,
                                                                                   'Uptake SD': sd,
                                                                                   'Uptake_SD_corr': sd_corr,
                                                                                   'Sequence': seq,
                                                                                   'Fragment':'C'}])
                                    else:
                                        if (srt_y[i] + 1) not in new_dict.keys():
                                            new_dict[(srt_y[i] + 1)] = {j: [{'Uptake': uptake,
                                                                             'Uptake SD': sd,
                                                                             'Sequence': seq,
                                                                             'Fragment': 'C'}]}

                                        else:
                                            if int(j) not in new_dict[srt_y[i] + 1].keys():
                                                new_dict[srt_y[i] + 1][j] = [{'Uptake': uptake,
                                                                              'Uptake SD': sd,
                                                                              'Sequence': seq,
                                                                              'Fragment': 'C'}]
                                            else:
                                                new_dict[srt_y[i] + 1][j].extend([{'Uptake': uptake,
                                                                                   'Uptake SD': sd,
                                                                                   'Sequence': seq,
                                                                                   'Fragment': 'C'}])

                    # Find matching C termini
                    # Make dictionary with N values for C keys from x dict
                    y_dict = {}
                    for x in x_dict.keys():
                        for y in x_dict[x].keys():
                            if y not in y_dict.keys():
                                y_dict[y] = {x: deepcopy(x_dict[x][y])}
                            elif x not in y_dict[y].keys():
                                y_dict[y][x] = deepcopy(x_dict[x][y])
                            else:
                                print("Error: Multiple matching Peptides")
                                break

                    # Iterate N values in C keys to create new peptides
                    for y in y_dict.keys():
                        srt_x = sorted(y_dict[y])
                        if len(srt_x) > 1:
                            for i in range(0, len(srt_x) - 1):
                                for j in srt_x[i + 1:]:
                                    uptake = round(float(y_dict[y][int(srt_x[i])]['Uptake']), res) - round(float(y_dict[y][j]['Uptake']), res)
                                    sd = sqrt(round(float(y_dict[y][int(srt_x[i])]['Uptake SD']), res) ** 2 + round(float(y_dict[y][j]['Uptake SD']), res) ** 2)
                                    seq = ''.join(y_dict[y][int(srt_x[i])]['Sequence'][:(j - srt_x[i]+1)])
                                    if self.corrected == 'Yes':
                                        uptake_corr = round(float(y_dict[y][int(srt_x[i])]['Uptake_corr']), res) - round(float(y_dict[y][j]['Uptake_corr']), res)
                                        sd_corr = sqrt(round(float(y_dict[y][int(srt_x[i])]['Uptake_SD_corr']), res) ** 2 + round(float(y_dict[y][j]['Uptake_SD_corr']), res) ** 2)

                                        if int(srt_x[i]) not in new_dict.keys():
                                            new_dict[int(srt_x[i])] = {(j): [{'Uptake': uptake,
                                                                                'Uptake_corr': uptake_corr,
                                                                                'Uptake SD': sd,
                                                                                'Uptake_SD_corr': sd_corr,
                                                                                'Sequence': seq,
                                                                                'Fragment':'N'}]}
                                        else:
                                            if int(j) not in new_dict[int(srt_x[i])].keys():
                                                new_dict[int(srt_x[i])][(j)] = [{'Uptake': uptake,
                                                                                     'Uptake_corr': uptake_corr,
                                                                                     'Uptake SD': sd,
                                                                                     'Uptake_SD_corr': sd_corr,
                                                                                     'Sequence': seq,
                                                                                     'Fragment':'N'}]
                                            else:
                                                new_dict[int(srt_x[i])][(j)].extend([{'Uptake': uptake,
                                                                                        'Uptake_corr': uptake_corr,
                                                                                        'Uptake SD': sd,
                                                                                        'Uptake_SD_corr': sd_corr,
                                                                                        'Sequence': seq,
                                                                                        'Fragment':'N'}])
                                    else:
                                        if int(srt_x[i]) not in new_dict.keys():
                                            new_dict[int(srt_x[i])] = {(j): [{'Uptake': uptake,
                                                                                'Uptake SD': sd,
                                                                                'Sequence': seq,
                                                                                'Fragment': 'N'}]}
                                        else:
                                            if int(j) not in new_dict[int(srt_x[i])].keys():
                                                new_dict[int(srt_x[i])][(j)] = [{'Uptake': uptake,
                                                                                   'Uptake SD': sd,
                                                                                   'Sequence': seq,
                                                                                   'Fragment': 'N'}]
                                            else:
                                                new_dict[int(srt_x[i])][(j)].extend([{'Uptake': uptake,
                                                                                        'Uptake SD': sd,
                                                                                        'Sequence': seq,
                                                                                        'Fragment': 'N'}])

                    # Calculate Max Uptake Values and export
                    for i in new_dict.keys():
                        for j in new_dict[i].keys():
                            for k in new_dict[i][j]:
                                self.pepid = self.pepid+1
                                up_count = 0
                                for l in list(k['Sequence']):
                                    if l != "P":
                                        up_count += 1
                                if self.corrected == 'Yes':
                                    if k['Fragment'] == 'N':
                                        maxuptake = up_count - 1
                                    else:
                                        maxuptake = up_count
                                    recombine_dict[protein][state][exposure].append(
                                        {"PepID": None, "Protein": protein, "Start": i, "End": j, "size": len(k['Sequence']),
                                         "Sequence": k['Sequence'], 'Modification': "", 'Fragment': k['Fragment'],
                                         'MaxUptake': maxuptake, 'MHP': "", "State": state, "Exposure": float(exposure),
                                         'Center': "",
                                         'Center SD': "",
                                         'Uptake': k['Uptake'],
                                         'Uptake_corr': k['Uptake_corr'],
                                         'Uptake SD': k['Uptake SD'],
                                         'Uptake_SD_corr': k['Uptake_SD_corr'],
                                         'RT': "", 'RT SD': "",
                                         'Percent_Uptake_corr': "",
                                         'ExposureCorrection': expf,
                                         'PeptideCorrection': pepf})
                                else:
                                    if k['Fragment'] == 'N':
                                        maxuptake = up_count - 1
                                    else:
                                        maxuptake = up_count
                                    recombine_dict[protein][state][exposure].append(
                                        {"PepID": None, "Protein": protein, "Start": i, "End": j, "size": len(k['Sequence']),
                                         "Sequence": k['Sequence'], 'Modification': "", 'Fragment': k['Fragment'],
                                         'MaxUptake': maxuptake, 'MHP': "", "State": state, "Exposure": float(exposure),
                                         'Center': "",
                                         'Center SD': "",
                                         'Uptake': k['Uptake'],
                                         'Uptake SD': k['Uptake SD'],
                                         'RT': "",
                                         'RT SD': ""})
        self.recombine_dict = recombine_dict

    # Add peptide generated from OPS to master_csv and data_nest
    def applySegmentation(self):
        correct = self.corrected
        for protein in self.proteins:
            for state in self.recombine_dict[protein].keys():
                for exposure in self.recombine_dict[protein][state].keys():
                    for row in self.recombine_dict[protein][state][exposure]:
                        if correct == 'Yes':
                            row['Uptake_corr'] = float(mean(row['Uptake_corr']))
                            row['Uptake_SD_corr'] = float(mean(row['Uptake_SD_corr']))
                            row['Uptake'] = float(mean(row['Uptake']))
                            row['Uptake SD'] = float(mean(row['Uptake SD']))
                        else:
                            row['Uptake'] = float(mean(row['Uptake']))
                            row['Uptake SD'] = float(mean(row['Uptake SD']))
                        self.master_csv.append(row) # Added so data table can update 20180909
                        self.data_nest[protein][state][float(exposure)].append(row)
        if self.filetype == 'CSV':
            self.assignIDs_CSV()
        else:
            self.assignIDs_DNX()
        if self.stats_nest is not None:
            self.stats_nest = None
            self.getSignificance()
        self.recombined = True

    # Perform retention time prediction
    def predictRT(self,seq):
        return achrom.calculate_RT(seq, self.seq_to_rt)

    # Export data to spreadsheet
    def exportData(self, export_name):
        if self.corrected == 'Yes':
            fieldnames = ['Protein', 'Start', 'End', 'PepID', 'size', 'Sequence', 'Modification', 'Fragment', 'MaxUptake',
                          'MHP', 'State', 'Exposure', 'Center', 'Center SD', 'Uptake', 'Uptake SD', 'RT', 'RT SD',
                          'Uptake_corr', 'Uptake_SD_corr', 'Percent_Uptake_corr', 'ExposureCorrection', 'PeptideCorrection']
        else:
            fieldnames = ['Protein', 'Start', 'End', 'PepID', 'size', 'Sequence', 'Modification', 'Fragment', 'MaxUptake',
                          'MHP', 'State', 'Exposure', 'Center', 'Center SD', 'Uptake', 'Uptake SD', 'RT', 'RT SD']
        with open(export_name, 'w') as name:
            writer = csv.DictWriter(name, fieldnames=fieldnames, extrasaction='ignore')
            writer.writeheader()
            for row in sorted(sorted(sorted(sorted(sorted(self.master_csv,
                                                          key=itemgetter('Exposure')),
                                                   key=itemgetter('State')),
                                            key=itemgetter('End')),
                                     key=itemgetter('Start')),
                              key=itemgetter('Protein')):
                writer.writerow(row)

    # Generate and Save PML script
    def exportPymol(self, corrected, state1, state2, protein, exposure, pyprot, pychain, type, color, range, export_name):
        if corrected == 'Yes':
            upval = "Uptake_corr"
            sdval = "Uptake_SD_corr"
        else:
            upval = "Uptake"
            sdval = "Uptake SD"
        if type in ['Uptake','SD','Fractional Uptake','Fractional SD']:
            if exposure == 'Ave':
                exposures = self.data_nest[protein][state1].keys()
                self.assignData(protein, state1, min(exposures))
            else:
                exposure = float(exposure)
                self.assignData(protein, state1, exposure)
        elif type in ['Fractional Uptake','Fractional Uptake Diff']:
            if exposure == 'Ave':
                exposures = list(
                    set(self.data_nest[protein][state1].keys()) & set(self.data_nest[protein][state2].keys()))
                self.assignData(protein, state1, min(exposures), state2)
            else:
                exposure = float(exposure)
                self.assignData(protein, state1, exposure, state2)
        else:
            print('Type Error')
        with open(export_name, 'w') as pmlfile:
            pmlfile.write('#Set a default B-Factor value for all relevant atoms\n')
            pmlfile.write('alter /' + pyprot + '//' + pychain + ', b=-99 \n')
            pmlfile.write('\n')
            pmlfile.write('#Set B-Factor for any residue atoms we have data for\n')
            values = []
            sizes = []
            if type == 'Uptake':
                for i in self.seq_nest:
                    if i['assigned'] == 1:
                        if exposure == 'Ave':
                            pepuptake = []
                            for e in sorted(exposures)[1:]:
                                pep = [n for n in self.data_nest[protein][state1][e] if
                                           n["Start"] == i['pepmin'] and n['End'] == i['pepmax']][0]
                                pepuptake.append(float(pep[upval]))
                            pepuptake = mean(pepuptake)
                        else:
                            pep = next(n for n in self.data_nest[protein][state1][exposure] if
                                       n["Start"] == i['pepmin'] and n['End'] == i['pepmax'])
                            pepuptake = pep[upval]
                        b = float(pepuptake)
                        values.append(b)
                        sizes.append(int(pep['size']))
                        pmlfile.write('alter /' + pyprot + '//' + pychain + '/' + str(i['pos']) + ', b=' + str(
                            round(b, 4)) + '\n')
                if range == 'Fit to Data':
                    minimum = min(values)
                    maximum = max(values)
                elif range == 'Full Range':
                    minimum = 0
                    maximum = max(sizes)
                elif range == 'Half Range':
                    minimum = 0
                    maximum = max(sizes) / 2
                else:
                    minimum = min(values)
                    maximum = max(values)
            elif type == 'SD':
                for i in self.seq_nest:
                    if i['assigned'] == 1:
                        if exposure == 'Ave':
                            pepsd = []
                            for e in sorted(exposures)[1:]:
                                pep = [n for n in self.data_nest[protein][state1][e] if
                                           n["Start"] == i['pepmin'] and n['End'] == i['pepmax']][0]
                                pepsd.append(float(pep[sdval]))
                            pepsd = mean(pepsd)
                        else:
                            pep = next(n for n in self.data_nest[protein][state1][exposure] if
                                       n["Start"] == i['pepmin'] and n['End'] == i['pepmax'])
                            pepsd = pep[sdval]
                        b = float(pepsd)
                        values.append(b)
                        pmlfile.write('alter /' + pyprot + '//' + pychain + '/' + str(i['pos']) + ', b=' + str(
                            round(b, 4)) + '\n')
                if range == 'Fit to Data':
                    minimum = min(values)
                    maximum = max(values)
                elif range == 'Full Range':
                    minimum = 0
                    maximum = max(values)
                elif range == 'Half Range':
                    minimum = 0
                    maximum = max(values)/2
                else:
                    minimum = min(values)
                    maximum = max(values)
            elif type == 'Uptake Diff':
                for i in self.seq_nest:
                    if i['assigned'] == 1:
                        if exposure == 'Ave':
                            b = []
                            for e in sorted(exposures)[1:]:
                                pep1 = [n[upval] for n in self.data_nest[protein][state1][e] if
                                            n["Start"] == i['pepmin'] and n['End'] == i['pepmax']][0]
                                pep2 = [n[upval] for n in self.data_nest[protein][state2][e] if
                                            n["Start"] == i['pepmin'] and n['End'] == i['pepmax']][0]
                                b.append((float(pep1) - float(pep2)))
                            b = mean(b)
                        else:
                            pep1 = next(n[upval] for n in self.data_nest[protein][state1][exposure] if
                                        n["Start"] == i['pepmin'] and n['End'] == i['pepmax'])
                            pep2 = next(n[upval] for n in self.data_nest[protein][state2][exposure] if
                                        n["Start"] == i['pepmin'] and n['End'] == i['pepmax'])
                            b = (float(pep1) - float(pep2))
                        values.append(b)
                        pmlfile.write('alter /' + pyprot + '//' + pychain + '/' + str(i['pos']) + ', b=' + str(
                            round(b, 4)) + '\n')
                if range == 'Fit to Data':
                    minimum = min(values)
                    maximum = max(values)
                elif range == 'Full Range':
                    minimum = -1
                    maximum = 1
                elif range == 'Half Range':
                    minimum = -0.5
                    maximum = 0.5
                else:
                    minimum = -1
                    maximum = 1
            elif type == 'Fractional Uptake':
                for i in self.seq_nest:
                    if i['assigned'] == 1:
                        if exposure == 'Ave':
                            pepuptake = []
                            for e in sorted(exposures)[1:]:
                                pep = [n for n in self.data_nest[protein][state1][e] if
                                           n["Start"] == i['pepmin'] and n['End'] == i['pepmax']][0]
                                pepuptake.append(float(pep[upval]))
                            pepuptake = mean(pepuptake)
                        else:
                            pep = next(n for n in self.data_nest[protein][state1][exposure] if
                                       n["Start"] == i['pepmin'] and n['End'] == i['pepmax'])
                            pepuptake = pep[upval]
                        maxuptake = pep["MaxUptake"]
                        b = float(pepuptake) / float(maxuptake)
                        values.append(b)
                        pmlfile.write('alter /' + pyprot + '//' + pychain + '/' + str(i['pos']) + ', b=' + str(
                            round(b, 4)) + '\n')
                if range == 'Fit to Data':
                    minimum = min(values)
                    maximum = max(values)
                elif range == 'Full Range':
                    minimum = 0
                    maximum = 1
                elif range == 'Half Range':
                    minimum = 0
                    maximum = 0.5
                else:
                    minimum = 0
                    maximum = 1
            elif type == 'Fractional SD (RSD)':
                for i in self.seq_nest:
                    if i['assigned'] == 1:
                        if exposure == 'Ave':
                            b = []
                            for e in sorted(exposures)[1:]:
                                pep = [n for n in self.data_nest[protein][state1][e] if
                                           n["Start"] == i['pepmin'] and n['End'] == i['pepmax']][0]
                                b.append(float(pep[sdval])*100 / float(pep[upval]))
                            b = mean(b)
                        else:
                            pep = next(n for n in self.data_nest[protein][state1][exposure] if
                                       n["Start"] == i['pepmin'] and n['End'] == i['pepmax'])
                            b = float(pep[sdval])*100 / float(pep[upval])
                        values.append(b)
                        pmlfile.write('alter /' + pyprot + '//' + pychain + '/' + str(i['pos']) + ', b=' + str(
                            round(b, 4)) + '\n')
                if range == 'Fit to Data':
                    minimum = min(values)
                    maximum = max(values)
                elif range == 'Full Range':
                    minimum = 0
                    maximum = max(values)
                elif range == 'Half Range':
                    minimum = 0
                    maximum = max(values)/2
                else:
                    minimum = min(values)
                    maximum = max(values)
            elif type == 'Fractional Uptake Diff':
                for i in self.seq_nest:
                    if i['assigned'] == 1:
                        if exposure == 'Ave':
                            b = []
                            for e in sorted(exposures)[1:]:
                                pep1 = [n[upval] for n in self.data_nest[protein][state1][e] if
                                            n["Start"] == i['pepmin'] and n['End'] == i['pepmax']][0]
                                pep2 = [n[upval] for n in self.data_nest[protein][state2][e] if
                                            n["Start"] == i['pepmin'] and n['End'] == i['pepmax']][0]
                                maxuptake = [n['MaxUptake'] for n in self.data_nest[protein][state1][exposure] if
                                                 n["Start"] == i['pepmin'] and n['End'] == i['pepmax']][0]
                                b.append((float(pep1) - float(pep2))/ float(maxuptake))
                            b = mean(b)
                        else:
                            pep1 = next(n[upval] for n in self.data_nest[protein][state1][exposure] if
                                        n["Start"] == i['pepmin'] and n['End'] == i['pepmax'])
                            maxuptake = next(n['MaxUptake'] for n in self.data_nest[protein][state1][exposure] if
                                             n["Start"] == i['pepmin'] and n['End'] == i['pepmax'])
                            pep2 = next(n[upval] for n in self.data_nest[protein][state2][exposure] if
                                        n["Start"] == i['pepmin'] and n['End'] == i['pepmax'])
                            b = (float(pep1) - float(pep2)) / float(maxuptake)
                        values.append(b)
                        pmlfile.write('alter /' + pyprot + '//' + pychain + '/' + str(i['pos']) + ', b=' + str(
                            round(b, 4)) + '\n')
                if range == 'Fit to Data':
                    minimum = min(values)
                    maximum = max(values)
                elif range == 'Full Range':
                    minimum = -1
                    maximum = 1
                elif range == 'Half Range':
                    minimum = -0.5
                    maximum = 0.5
                else:
                    minimum = -1
                    maximum = 1
            else:
                print('Type Error')
            pmlfile.write('select missing, /' + pyprot + '//' + str(pychain) + ' and b = -99\n')
            pmlfile.write('spectrum b, ' + str(color) + ', minimum=' + str(minimum) + ', maximum=' + str(
                maximum) + ', selection=/'
                          + pyprot + '//' + str(pychain) + '\n')
            pmlfile.write('color black, missing\n')

    # Export Outlier Statistics to spreadsheet
    def exportStats(self, export_name):
        fieldnames = ['Protein', 'Start', 'End', 'Sequence', 'Modification', 'Fragment', 'MaxUptake', 'MHP', 'State',
                      'Exposure', 'Center', 'CenterSD', 'RawID', 'z', 'MZ', 'PPM', 'PPM_ave', 'PPM_sd', 'PPM_5',
                      'PPM_10', 'PPM_25', 'RT', 'RT_ave', 'RT_sd', 'RT_zscore', 'Mob', 'Mob_ave', 'Mob_sd',
                      'Mob_zscore']
        cluster = self.getDnxCluster([])
        with open(export_name, 'w') as name:
            writer = csv.DictWriter(name, fieldnames=fieldnames, extrasaction='ignore')
            writer.writeheader()
            sortedcluster = copy(sorted(
                sorted(sorted(sorted(sorted(sorted(sorted(cluster, key=itemgetter('z')),
                                                   key=itemgetter('RawID')),
                                            key=itemgetter('Exposure')),
                                     key=itemgetter('State')),
                              key=itemgetter('Length')),
                       key=itemgetter('Start')),
                key=itemgetter('Protein')))
            for row in sortedcluster:
                stats = self.file_nest[row['PepID']]['States'][row['State']][row['Exposure']][row['RawID']][row['z']]
                newrow = {}
                newrow.update({'Protein':row['Protein'],
                               'Start':int(row['Start']) + 1,
                               'End':int(row['Start']) + int(row['Length']),
                               'Sequence':row['Sequence'],
                               'Modification':row['Modification'],
                               'Fragment':row['Fragment'],
                               'MaxUptake':row['MaxUptake'],
                               'MHP':row['MHP'],
                               'State':row['State'],
                               'Exposure':row['Exposure'],
                               'Center':row['Center'],
                               'CenterSD':row['CenterSD'],
                               'RawID':row['RawID'],
                               'z':row['z'],
                               'MZ':stats['data'],
                               'PPM':stats['ppm']['data'],
                               'PPM_ave':mean(stats['ppm']['data']),
                               'PPM_sd':std(stats['ppm']['data']),
                               'RT':stats['rt']['data'],
                               'RT_ave':mean(stats['rt']['data']),
                               'RT_sd':std(stats['rt']['data']),
                               'Mob':stats['drift']['data'],
                               'Mob_ave':mean(stats['drift']['data']),
                               'Mob_sd':std(stats['drift']['data'])})
                writer.writerow(newrow)

    # Export Significance Statistics to spreadsheet
    def exportSigs(self, export_name):
        newfields = [str(i[0]) + '_' + str(i[1]) for i in list(combinations(self.states, 2))]

        fieldnames = ['Protein', 'Start', 'End', 'size', 'Sequence', 'Modification', 'Fragment', 'MaxUptake',
                      'MHP', 'State', 'Exposure', 'Center', 'Center SD', 'Uptake', 'Uptake SD', 'RT', 'RT SD',
                      'LeveneF', 'LeveneP', 'AnovaF', 'AnovaP']
        fieldnames.extend(newfields)
        with open(export_name, 'w') as name:
            writer = csv.DictWriter(name, fieldnames=fieldnames, extrasaction='ignore')
            writer.writeheader()
            for row in self.master_csv:
                sigs = self.stats_nest[row['PepID']][row['Exposure']]['Stats']
                row.update({'LeveneF': sigs['Levene']['F'],
                            'LeveneP': sigs['Levene']['p'],
                            'AnovaF': sigs['Anova']['F'],
                            'AnovaP':sigs['Anova']['p']})
                for i in sigs['Tukey']:
                    name1 = i['group1'] + '_' + i['group2']
                    name2 = i['group2'] + '_' + i['group1']
                    if name1 in fieldnames:
                        row[name1] = i['p-value']
                    elif name2 in fieldnames:
                        row[name2] = i['p-value']
                writer.writerow(row)

# Private Function from statsmodels package, to be able to calculate CIs from summary statistics
def simultaneous_ci(q_crit, var, groupnobs): # from statmodels package
    """Compute simultaneous confidence intervals for comparison of means.

    q_crit value is generated from tukey hsd test. Variance is considered
    across all groups. Returned halfwidths can be thought of as uncertainty
    intervals around each group mean. They allow for simultaneous
    comparison of pairwise significance among any pairs (by checking for
    overlap)

    Parameters
    ----------
    q_crit : float
        The Q critical value studentized range statistic from Tukey's HSD
    var : float
        The group variance
    groupnobs : array-like object
        Number of observations contained in each group.
    pairindices : tuple of lists, optional
        Indices corresponding to the upper triangle of matrix. Computed
        here if not supplied

    Returns
    -------
    halfwidths : ndarray
        Half the width of each confidence interval for each group given in
        groupnobs

    See Also
    --------
    MultiComparison : statistics class providing significance tests
    tukeyhsd : among other things, computes q_crit value

    References
    ----------
    .. [*] Hochberg, Y., and A. C. Tamhane. Multiple Comparison Procedures.
           Hoboken, NJ: John Wiley & Sons, 1987.)
    """
    # Set initial variables
    ng = len(groupnobs)

    # Compute dij for all pairwise comparisons ala hochberg p. 95
    gvar = var / groupnobs
    pairindices = triu_indices(ng, 1)
    d12 = sqrt(gvar[pairindices[0]] + gvar[pairindices[1]])

    # Create the full d matrix given all known dij vals
    d = zeros((ng, ng))
    d[pairindices] = d12
    d = d + d.conj().T

    # Compute the two global sums from hochberg eq 3.32
    sum1 = sum(d12)
    sum2 = sum(d, axis=0)

    if (ng > 2):
        w = ((ng-1.) * sum2 - sum1) / ((ng - 1.) * (ng - 2.))
    else:
        w = sum1 * ones((2, 1)) / 2.

    return (q_crit / sqrt(2))*w