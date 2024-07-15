import numpy as np
import copy
import json


def Sort_LG(_Matrix_LG):
    Out_Matrix_LG = []
    for j in range(0, len(_Matrix_LG)):
        _Matrix = _Matrix_LG[j]
        Complexes = {}
        for i in range(0, len(_Matrix)):
            if _Matrix[i] not in Complexes:
                Complexes[_Matrix[i]] = [i]
            else:
                Complexes[_Matrix[i]].append(i)
        Temp_LG_List = ['a' for i in range(0, len(_Matrix))]
        To_Count = 0
        To_fill = 0
        Done_complexes = {}
        while len(Done_complexes) < len(Complexes):
            for keys in Complexes:
                if keys not in Done_complexes:
                    if To_Count in Complexes[keys]:
                        for each_index in Complexes[keys]:
                            Temp_LG_List[each_index] = To_fill
                        To_fill = To_fill + 1
                        Done_complexes[keys] = ''
                    else:
                        pass
                else:
                    pass
            To_Count = To_Count + 1
        Out_Matrix_LG.append(Temp_LG_List)
    return Out_Matrix_LG

def Matrix2String01_LG(_Matrix):
    '''Convert numpy array into 1-D string representation of Logic Gate Matrix'''
    _Matrix = np.array(_Matrix)
    OutString01 = ''
    for i in range(0, _Matrix.shape[0]):
        for j in range(0, _Matrix.shape[1]):
            OutString01 = OutString01 + str(_Matrix[i][j])
    return OutString01

def Matrix2String01_LG_Expanded(_Matrix):
    '''Convert numpy array into 1-D string representation of Logic Gate Matrix'''
    _Matrix = np.array(_Matrix)
    OutString01 = ''
    for i in range(0, _Matrix.shape[0]):
        for j in range(0, _Matrix.shape[1]):
            OutString01 = OutString01 + str(_Matrix[i][j]) + ','
    OutString01 = OutString01[:-1]
    return OutString01

def LogicGatesString2Matrix(string):
    '''Convert 1-D string representation of Logic Gate Matrix into numpy array for fast matrix manipulations'''
    outmatrix = np.random.randint(3, 4, (int(len(string)/2), 2))
    stringindex = 0
    for i in range(0, outmatrix.shape[0]):
        for j in range(0, outmatrix.shape[1]):
            outmatrix[i][j] = int(string[stringindex])
            stringindex = stringindex + 1
    return outmatrix

def LogicGatesString2Matrix_Expanded(string):
    '''Convert 1-D string representation of Logic Gate Matrix into numpy array for fast matrix manipulations'''
    string_matrix = string.split(',')
    SQRTLEN = int(np.sqrt(len(string_matrix)))
    outmatrix = np.random.randint(8, 9, (SQRTLEN, SQRTLEN))
    stringindex = 0
    for i in range(0, outmatrix.shape[0]):
        for j in range(0, outmatrix.shape[1]):
            #print(stringindex)
            outmatrix[i][j] = int(string_matrix[stringindex])
            stringindex = stringindex + 1
    return outmatrix

def GetCorespondingMatrix(List012):
    '''Convert 1-D array representation of Weighted Adjacency Matrix into numpy array for fast matrix manipulations'''
    SQRTLEN = int(np.sqrt(len(List012)))
    OutMatrix = np.random.randint(8, 9, (2, SQRTLEN, SQRTLEN))
    for i in range(0, len(List012)):
        if List012[i] == 0:
            OutMatrix[0][i//SQRTLEN][i % SQRTLEN] = 0
            OutMatrix[1][i//SQRTLEN][i % SQRTLEN] = 0
        elif List012[i] == 1:
            OutMatrix[0][i//SQRTLEN][i % SQRTLEN] = 1
            OutMatrix[1][i//SQRTLEN][i % SQRTLEN] = 0
        elif List012[i] == 2:
            OutMatrix[0][i//SQRTLEN][i % SQRTLEN] = 0
            OutMatrix[1][i//SQRTLEN][i % SQRTLEN] = 1
        else:
            raise Exception('List012 must just have 012!\n{}'.format(List012))
    return OutMatrix


def String012ToMatrix(List012):
    '''Convert 1-D string representation of Weighted Adjacency Matrix into numpy array for fast matrix manipulations'''
    # String '012' to Matrix
    SQRTLEN = int(np.sqrt(len(List012)))
    OutMatrix = np.random.randint(8, 9, (2, SQRTLEN, SQRTLEN))
    for i in range(0, len(List012)):
        if List012[i] == '0':
            OutMatrix[0][i//SQRTLEN][i % SQRTLEN] = 0
            OutMatrix[1][i//SQRTLEN][i % SQRTLEN] = 0
        elif List012[i] == '1':
            OutMatrix[0][i//SQRTLEN][i % SQRTLEN] = 1
            OutMatrix[1][i//SQRTLEN][i % SQRTLEN] = 0
        elif List012[i] == '2':
            OutMatrix[0][i//SQRTLEN][i % SQRTLEN] = 0
            OutMatrix[1][i//SQRTLEN][i % SQRTLEN] = 1
        else:
            raise Exception('List012 must just have 012!\n{}'.format(List012))
    return OutMatrix


def GetCorespondingMatrix_LG(List01):
    '''Convert 1-D array representation of Logic Gate Matrix into numpy array for fast matrix manipulations'''
    SQRTLEN = int(len(List01)/2)
    OutMatrix = np.random.randint(8, 9, (SQRTLEN, 2))
    for i in range(0, len(List01)):
        OutMatrix[i//2][i % 2] = int(List01[i])
    return OutMatrix


def ConfigurationTo012(ConfigurationMatrix):
    '''Convert configuration matrix to string'''
    OutString = ''
    for j in range(0, ConfigurationMatrix.shape[1]):
        for z in range(0, ConfigurationMatrix.shape[2]):
            if ConfigurationMatrix[0][j][z] == 0 and ConfigurationMatrix[1][j][z] == 0:
                OutString = OutString + '0'
            elif ConfigurationMatrix[0][j][z] == 1 and ConfigurationMatrix[1][j][z] == 0:
                OutString = OutString + '1'
            elif ConfigurationMatrix[0][j][z] == 0 and ConfigurationMatrix[1][j][z] == 1:
                OutString = OutString + '2'
            else:
                raise Exception('Configuration (1,1)!')
    return OutString


def maskmean(array_):
    '''mean of array_ with entries called 'mask' removed'''
    out_mean = []
    for each_array in array_.T:
        temp_mean = []
        for each in each_array[0]:
            if each == 'mask':
                pass
            else:
                temp_mean.append(float(each))
        out_mean.append(np.nanmean(temp_mean))
    return out_mean

def GetmRNASearchSpace(TotalGeneNum, RunNum, BoundList):
    '''obtain the coordinates of evenly-distributed vertices in a n-dimensional space'''
    BoundDeltaList = []
    for i in range(0, len(BoundList)):
        BoundDeltaList.append(abs(BoundList[i][0] - BoundList[i][1]))
    if len(BoundDeltaList) == TotalGeneNum:
        LianJi = 1.0
        for i in range(0, TotalGeneNum):
            LianJi = LianJi * BoundDeltaList[i]
        Delta = round((LianJi / RunNum) ** (1 / TotalGeneNum))
        if Delta == 0:
            Delta = 1
        else:
            pass
    else:
        raise Exception('len(BoundDeltaList) != TotalGeneNum')
    Delta = int(Delta)

    ExpandedList = []
    OutList = []
    for i in range(0, len(BoundList)):
        ExpandedList.append([x for x in range(min(BoundList[i]), max(BoundList[i]) + 1, Delta)])
    CounterList = [0 for i in range(0, len(ExpandedList))]
    CounterListEndState = [max(ExpandedList[p]) for p in range(0, len(ExpandedList))]
    while True:
        OutList.append(copy.deepcopy(CounterList))
        if CounterList == CounterListEndState:
            break
        else:
            pass
        for i in range(0, len(CounterList)):
            if CounterList[i] >= CounterListEndState[i]:
                CounterList[i] = 0
            else:
                CounterList[i] = CounterList[i] + Delta
                break
    return np.array(OutList)

def json2ea(json_data):
    '''Convert .json to input format for EA'''
    # Parse JSON data
    data = json.loads(json_data)

    # Extract nodes and edges from JSON
    nodes = {node["id"]: node["f0"] for node in data["nodes"]}
    edges = [(edge["source"], edge["target"], edge["label"], edge["style"]) for edge in data["edges"]]

    # Identify unique node IDs
    Order_of_genes = list(nodes.keys())
    num_nodes = len(Order_of_genes)

    # Initialize adjacency matrix
    adjacency_matrix = [[0] * num_nodes for _ in range(num_nodes)]

    # Get the adjacency matrix
    for edge in edges:
        source_index = Order_of_genes.index(edge[0])
        target_index = Order_of_genes.index(edge[1])
        if edge[3][0] == 'triangle':
            adjacency_matrix[source_index][target_index] = 1
        elif edge[3][0] == 'tee':
            adjacency_matrix[source_index][target_index] = 2
        else:
            adjacency_matrix[source_index][target_index] = 0
    AM = ''
    for each_x in adjacency_matrix:
        for each_y in each_x:
            AM = AM + str(each_y)

    # Get the f0
    f0 = []
    for node in Order_of_genes:
        f0.append(nodes[node])
    f0 = str(f0)

    temp_LG = {}

    for edge in edges:
        for each in Order_of_genes:
            if edge[1] == each and edge[3][1] == 'dashed':
                if edge[1] not in temp_LG:
                    temp_LG[edge[1]] = {edge[0]: edge[2]}
                else:
                    temp_LG[edge[1]].update({edge[0]: edge[2]})

    LG = []
    for each in Order_of_genes:
        if each not in temp_LG:
            LG.append([i for i in range(0, num_nodes)])
        else:
            LG.append(['a' for i in range(0, num_nodes)])
            counter = 0
            copy_temp_LG = copy.deepcopy(temp_LG[each])
            for gene in Order_of_genes:
                if gene in copy_temp_LG:
                    keys_to_remove = []
                    for keys in temp_LG[each]:
                        if (gene in temp_LG[each]) and (temp_LG[each][keys] == temp_LG[each][gene]):
                            LG[-1][Order_of_genes.index(keys)] = counter
                            keys_to_remove.append(keys)
                        else:
                            pass
                    for key_to_remove in keys_to_remove:
                        del temp_LG[each][key_to_remove]
                    if keys_to_remove == []:
                        pass
                    else:
                        counter = counter + 1
                else:
                    LG[-1][Order_of_genes.index(gene)] = counter
                    counter = counter + 1

    LG_ = ''
    for each_i in LG:
        for each_j in each_i:
            LG_ = LG_ + str(each_j) + ','
    LG_ = LG_[:-1]

    return AM, f0, LG_, '\t'.join(Order_of_genes)

def ea2json(AM, LG, f0, Order_of_genes, Vmax=''):
    '''Convert EA output to .json'''
    if Vmax == '':
        Vmax = [1 for i in range(0, len(Order_of_genes))]
    else:
        Vmax = eval(Vmax)

    json = ''
    # Add the nodes
    f0 = eval(f0)

    LG = LG.split(',')

    json = json + '{\n  "nodes": [\n'
    for i in range(0, len(Order_of_genes)):
        json = json + '\t{\n\t  ' + '"id": "{}",\n'.format(Order_of_genes[i]) + '\t  "label": "{}",\n'.format(Order_of_genes[i]) + '\t  "sua7Occupancy": {},\n'.format(Vmax[i]) + '\t  "f0": {}'.format(f0[i])+ '\n  \t},\n'
    json = json[:-2]

    # Add the edges
    json = json + '\n  ],\n  "edges": [\n'
    num_nodes = len(f0)
    for i in range(0, len(AM)):
        if AM[i] == '0':
            continue
        elif AM[i] == '1':
            LG_slide = LG[num_nodes*(i%num_nodes):num_nodes*(1+i%num_nodes)]
            #print(i, LG_slide, i//num_nodes)
            if LG_slide.count(LG_slide[i//num_nodes]) > 1:
                json = json + '\t{\n\t  ' + '"source": "{}",\n'.format(Order_of_genes[i//num_nodes]) + '\t  "target": "{}",\n'.format(Order_of_genes[i%num_nodes]) + '\t  "label": "{}",\n'.format(LG_slide[i//num_nodes]) + '\t  "style": [\n\t\t"dashed",\n\t\t"triangle"\n\t  ]'+ '\n  \t},\n'
            else:
                json = json + '\t{\n\t  ' + '"source": "{}",\n'.format(Order_of_genes[i//num_nodes]) + '\t  "target": "{}",\n'.format(Order_of_genes[i%num_nodes]) + '\t  "label": "",\n' + '\t  "style": [\n\t\t"solid",\n\t\t"triangle"\n\t  ]'+ '\n  \t},\n'
        else:
            LG_slide = LG[num_nodes*(i%num_nodes):num_nodes*(1+i%num_nodes)]
            #print(i, LG_slide, i//num_nodes)
            if LG_slide.count(LG_slide[i//num_nodes]) > 1:
                json = json + '\t{\n\t  ' + '"source": "{}",\n'.format(Order_of_genes[i//num_nodes]) + '\t  "target": "{}",\n'.format(Order_of_genes[i%num_nodes]) + '\t  "label": "{}",\n'.format(LG_slide[i//num_nodes]) + '\t  "style": [\n\t\t"dashed",\n\t\t"tee"\n\t  ]'+ '\n  \t},\n'
            else:
                json = json + '\t{\n\t  ' + '"source": "{}",\n'.format(Order_of_genes[i//num_nodes]) + '\t  "target": "{}",\n'.format(Order_of_genes[i%num_nodes]) + '\t  "label": "",\n' + '\t  "style": [\n\t\t"solid",\n\t\t"tee"\n\t  ]'+ '\n  \t},\n'
    json = json[:-2]
    json = json + '\n  ]\n}'
    return json
