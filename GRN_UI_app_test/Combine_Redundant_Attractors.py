import numpy as np
from distance_functions import *
import sys

def to_combine(_matrix, AttractorDistance_Threshold, Identical_Expression_Threshold):
    '''Determine if a set of transcriptional profiles have two whose attractor distance is below threshold.'''
    '''Also return two transcriptional profiles if there's one and only one variable gene.'''
    _matrix = np.array(_matrix)
    for i in range(0, _matrix.shape[0]):
        for j in range(0, _matrix.shape[0]):
            if i == j:
                continue
            else:
                if GetAttractorDistance(_matrix[i], _matrix[j], np.max(_matrix, axis=0)) <= AttractorDistance_Threshold:
                    return False, min(i, j), max(i, j)
                else:
                    pass

                temp_vector = _matrix[j] - _matrix[i]
                for k in (0, len(temp_vector)):
                    if k == len(temp_vector):
                        if GetAttractorDistance(temp_vector[:-1], np.array([0 for repeat in range(0, len(temp_vector)-1)]), np.max(_matrix, axis=0)[:-1]) <= Identical_Expression_Threshold:
                            return False, min(i, j), max(i, j)
                        else:
                            pass
                    elif k == 0:
                        if GetAttractorDistance(temp_vector[1:], np.array([0 for repeat in range(0, len(temp_vector)-1)]), np.max(_matrix, axis=0)[1:]) <= Identical_Expression_Threshold:
                            return False, min(i, j), max(i, j)
                        else:
                            pass
                    else:
                        if GetAttractorDistance(temp_vector[:k]+temp_vector[k+1:], np.array([0 for repeat in range(0, len(temp_vector)-1)]), np.max(_matrix, axis=0)[:k]+np.max(_matrix, axis=0)[k+1:]) <= Identical_Expression_Threshold:
                            return False, min(i, j), max(i, j)
                        else:
                            pass

    return True, min(i, j), max(i, j)

def Combine_Redundant_Attractors(file_or_matrix=0, Matrix=[], Path_input='', AttractorDistance_Threshold=0.1, Identical_Expression_Threshold=0.05):
    '''Combine redundant attractors according to criterion specified in to_combine.'''

    TranscriptionPofiles = []

    '''Specify a genotype for combining.'''
    if file_or_matrix == 1:
        RNAseq_data = open(Path_input, 'r')
        for line in RNAseq_data:
            line_split = line.split()
            if line_split[0] != '-1':
                continue
            else:
                for i in range(0, len(line_split)):
                    line_split[i] = 100*float(line_split[i])
            TranscriptionPofiles.append(line_split[1:])
    else:
        TranscriptionPofiles = Matrix

    while True:
        _end, i, j = to_combine(TranscriptionPofiles, AttractorDistance_Threshold, Identical_Expression_Threshold)
        print(str(len(TranscriptionPofiles)))
        sys.stdout.flush() 
        if _end:
            break
        else:
            new_profile = (np.array(TranscriptionPofiles[i]) + np.array(TranscriptionPofiles[j])) * 0.5
            del TranscriptionPofiles[j]
            del TranscriptionPofiles[i]
            TranscriptionPofiles.insert(i, list(new_profile))
    return TranscriptionPofiles

