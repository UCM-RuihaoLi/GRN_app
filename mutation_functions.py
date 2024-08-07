from utility_functions import *
from distance_functions import *
import numpy as np
import random


def HammingMutation_LG(_HammingDistance, _Matrix, _globalmutationdirection, _samplesize):
    '''Hamming mutation function for Logic Gate Matrix'''
    Temp_String = []
    _S = Matrix2String01_LG(_Matrix)
    for i in range(0, len(_S)):
        Temp_String.append(int(_S[i]))
    MutationIndexList = [i for i in range(0, len(Temp_String))]
    Counteri = 0
    while Counteri < _HammingDistance:
        Index_to_Mutate = random.choice(MutationIndexList)
        FrequencyList = [_samplesize-_globalmutationdirection[str(Index_to_Mutate)][0],
                         _samplesize-_globalmutationdirection[str(Index_to_Mutate)][1]]  # Get [30,30,40]
        ProbList = []
        for ai in range(0, len(FrequencyList)):
            for aj in range(0, FrequencyList[ai]):
                if ai == 0:
                    ProbList.append(0)
                elif ai == 1:
                    ProbList.append(1)
                else:
                    raise Exception('lenFrequencyList is Wrong!')
        ProbList = np.array(ProbList)
        ProbList = ProbList.flatten()
        MutationResult = np.random.choice(ProbList)
        if Temp_String[Index_to_Mutate] == MutationResult:
            pass
        else:
            Temp_String[Index_to_Mutate] = MutationResult
            Counteri = Counteri + 1
            del MutationIndexList[MutationIndexList.index(Index_to_Mutate)]
    return GetCorespondingMatrix_LG(Temp_String)

def Mutation_LG(_HammingDistance, _Matrix_LG, _AR_List, PPI_='', prob=0.9):
    # Mutation for the expanded combinatorial LG.
    Counteri = 0
    if type(PPI_) == str:
        while Counteri < _HammingDistance:
            Gene_to_mutate = np.random.choice([i for i in range(0, len(_Matrix_LG[0]))])
            TF_to_mutate = np.random.choice([i for i in range(0, len(_Matrix_LG[0]))])
            if TF_to_mutate in _AR_List[0][Gene_to_mutate]:
                Acceptable_complex_to_mutate = []
                for each_TF in range(0, len(_Matrix_LG[0])):
                    if each_TF not in _AR_List[1][Gene_to_mutate] and each_TF != _Matrix_LG[Gene_to_mutate][TF_to_mutate]:
                        Acceptable_complex_to_mutate.append(each_TF)
                    else:
                        pass
            else:
                Acceptable_complex_to_mutate = []
                for each_TF in range(0, len(_Matrix_LG[0])):
                    if each_TF not in _AR_List[0][Gene_to_mutate] and each_TF != _Matrix_LG[Gene_to_mutate][TF_to_mutate]:
                        Acceptable_complex_to_mutate.append(each_TF)
                    else:
                        pass
            _Matrix_LG[Gene_to_mutate][TF_to_mutate] = np.random.choice(Acceptable_complex_to_mutate)
            Counteri = Counteri + 1
    else:
        while Counteri < _HammingDistance:
            Gene_to_mutate = np.random.choice([i for i in range(0, len(_Matrix_LG[0]))])
            TF_to_mutate = np.random.choice([i for i in range(0, len(_Matrix_LG[0]))])
            if TF_to_mutate in _AR_List[0][Gene_to_mutate]:
                Acceptable_complex_to_mutate = []
                for each_TF in range(0, len(_Matrix_LG[0])):
                    if each_TF not in _AR_List[1][Gene_to_mutate] and each_TF != _Matrix_LG[Gene_to_mutate][TF_to_mutate]:
                        TFs_in_complex = [index_ for index_, value_ in enumerate(_Matrix_LG[Gene_to_mutate]) if value_ == _Matrix_LG[Gene_to_mutate][TF_to_mutate]] # get the rest of the TFs in this complex
                        if any(PPI_[TF_to_mutate][ii_] == 1 for ii_ in TFs_in_complex): # the TF_to_mutate has PPI with any of the TFs in the complex
                            Acceptable_complex_to_mutate.extend([each_TF for _ in range(int(10*prob))])
                        else:
                            Acceptable_complex_to_mutate.append(each_TF)
                    else:
                        pass
            else:
                Acceptable_complex_to_mutate = []
                for each_TF in range(0, len(_Matrix_LG[0])):
                    if each_TF not in _AR_List[0][Gene_to_mutate] and each_TF != _Matrix_LG[Gene_to_mutate][TF_to_mutate]:
                        TFs_in_complex = [index_ for index_, value_ in enumerate(_Matrix_LG[Gene_to_mutate]) if value_ == _Matrix_LG[Gene_to_mutate][TF_to_mutate]] # get the rest of the TFs in this complex
                        if any(PPI_[TF_to_mutate][ii_] == 1 for ii_ in TFs_in_complex):  # the TF_to_mutate has PPI with any of the TFs in the complex
                            Acceptable_complex_to_mutate.extend([each_TF for _ in range(int(10 * prob))])
                        else:
                            Acceptable_complex_to_mutate.append(each_TF)
                    else:
                        pass
            _Matrix_LG[Gene_to_mutate][TF_to_mutate] = np.random.choice(Acceptable_complex_to_mutate)
            Counteri = Counteri + 1

    return Sort_LG(_Matrix_LG)


def HammingMutation(sys_input_ChIP,
                    _HammingDistance, _Matrix,
                    _globalmutationdirection, _samplesize, _Probability):
    '''Hamming mutation function for Adjacency Matrix'''
    # The shape of matrix is (2,x,x)
    Temp_String = []
    for i in range(0, _Matrix.shape[1]):
        for j in range(0, _Matrix.shape[2]):
            if _Matrix[0][i][j] == 0 and _Matrix[1][i][j] == 0:
                Temp_String.append(0)
            elif _Matrix[0][i][j] == 1 and _Matrix[1][i][j] == 0:
                Temp_String.append(1)
            elif _Matrix[0][i][j] == 0 and _Matrix[1][i][j] == 1:
                Temp_String.append(2)
            else:
                raise Exception('Matrix has something wrong!')
    if type(sys_input_ChIP) == str:
        MutationIndexList = [i for i in range(0, len(Temp_String))]
    else:
        MutationIndexList = []
        while len(MutationIndexList) <= 2:
            for ChIP_i in range(0, sys_input_ChIP.shape[0]):
                for ChIP_j in range(0, sys_input_ChIP.shape[1]):
                    if sys_input_ChIP[ChIP_i][ChIP_j] == 0:
                        if _Matrix[0][ChIP_i][ChIP_j] == 0 and _Matrix[1][ChIP_i][ChIP_j] == 0:
                            if np.random.uniform() < _Probability:
                                pass
                            else:
                                MutationIndexList.append(
                                    ChIP_i*sys_input_ChIP.shape[1]+ChIP_j)
                        else:
                            if np.random.uniform() >= _Probability:
                                pass
                            else:
                                MutationIndexList.append(
                                    ChIP_i*sys_input_ChIP.shape[1]+ChIP_j)
                    elif sys_input_ChIP[ChIP_i][ChIP_j] == 1:
                        if _Matrix[0][ChIP_i][ChIP_j] == 0 and _Matrix[1][ChIP_i][ChIP_j] == 0:
                            if np.random.uniform() >= _Probability:
                                pass
                            else:
                                MutationIndexList.append(
                                    ChIP_i*sys_input_ChIP.shape[1]+ChIP_j)
                        else:
                            if np.random.uniform() < _Probability:
                                pass
                            else:
                                MutationIndexList.append(
                                    ChIP_i*sys_input_ChIP.shape[1]+ChIP_j)
                    else:
                        raise Exception('ChIP data matrix is wrong!')
    Counteri = 0
    while Counteri < _HammingDistance:
        Index_to_Mutate = random.choice(MutationIndexList)
        del MutationIndexList[MutationIndexList.index(Index_to_Mutate)]
        if Temp_String[Index_to_Mutate] == 0:
            FrequencyList = [_samplesize-_globalmutationdirection[str(Index_to_Mutate)][1],
                             _samplesize-_globalmutationdirection[str(Index_to_Mutate)][2]]  # Get [30,30,40]
            ProbList = []
            for ai in range(0, len(FrequencyList)):
                for aj in range(0, FrequencyList[ai]):
                    if ai == 0:
                        ProbList.append(1)
                    elif ai == 1:
                        ProbList.append(2)
                    else:
                        raise Exception('lenFrequencyList is Wrong!')
            ProbList = np.array(ProbList)
            ProbList = ProbList.flatten()
            MutationResult = np.random.choice(ProbList)
            Temp_String[Index_to_Mutate] = MutationResult
        elif Temp_String[Index_to_Mutate] == 1 and Counteri < (_HammingDistance-1):
            FrequencyList = [_samplesize-_globalmutationdirection[str(Index_to_Mutate)][0],
                             _samplesize-_globalmutationdirection[str(Index_to_Mutate)][2]]  # Get [30,30,40]
            ProbList = []
            for ai in range(0, len(FrequencyList)):
                for aj in range(0, FrequencyList[ai]):
                    if ai == 0:
                        ProbList.append(0)
                    elif ai == 1:
                        ProbList.append(2)
                    else:
                        raise Exception('lenFrequencyList is Wrong!')
            ProbList = np.array(ProbList)
            ProbList = ProbList.flatten()
            MutationResult = np.random.choice(ProbList)
            Temp_String[Index_to_Mutate] = MutationResult
            if MutationResult == 2:
                Counteri = Counteri + 1
            else:
                pass
        elif Temp_String[Index_to_Mutate] == 2 and Counteri < (_HammingDistance-1):
            FrequencyList = [_samplesize-_globalmutationdirection[str(Index_to_Mutate)][0],
                             _samplesize-_globalmutationdirection[str(Index_to_Mutate)][1]]  # Get [30,30,40]
            ProbList = []
            for ai in range(0, len(FrequencyList)):
                for aj in range(0, FrequencyList[ai]):
                    if ai == 0:
                        ProbList.append(0)
                    elif ai == 1:
                        ProbList.append(1)
                    else:
                        raise Exception('lenFrequencyList is Wrong!')
            ProbList = np.array(ProbList)
            ProbList = ProbList.flatten()
            MutationResult = np.random.choice(ProbList)
            Temp_String[Index_to_Mutate] = MutationResult
            if MutationResult == 1:
                Counteri = Counteri + 1
            else:
                pass
        elif Counteri == (_HammingDistance-1) and Temp_String[Index_to_Mutate] != 0:
            Temp_String[Index_to_Mutate] = 0
        else:
            raise Exception('Temp_String has something wrong!')
        Counteri = Counteri + 1
    return GetCorespondingMatrix(Temp_String)


def Recombination(MatrixA, MatrixB):
    '''Exchange a slice between the two matrice given'''
    if MatrixA.shape == MatrixB.shape:
        Slice_X = sorted(np.random.randint(MatrixA.shape[1], size=(1, 2))[0])
        Slice_Y = sorted(np.random.randint(MatrixA.shape[2], size=(1, 2))[0])
        Slice_W = np.random.randint(2)
        # left closed right open
        A_Slice = MatrixA[Slice_W][Slice_X[0]: Slice_X[1] +
                                   1, Slice_Y[0]: Slice_Y[1]+1].copy()
        B_Slice = MatrixB[Slice_W][Slice_X[0]: Slice_X[1] +
                                   1, Slice_Y[0]: Slice_Y[1]+1].copy()
        MatrixA[Slice_W][Slice_X[0]: Slice_X[1] +
                         1, Slice_Y[0]: Slice_Y[1]+1] = B_Slice
        MatrixB[Slice_W][Slice_X[0]: Slice_X[1] +
                         1, Slice_Y[0]: Slice_Y[1]+1] = A_Slice
    else:
        raise Exception('The shape of MatrixA and MatrixB do not match!')
    return
