import copy
import numpy as np
import math
import random

from distance_functions import *
from mutation_functions import *
from dynamics import *
# non-overlapping generations in a haploid infinite population, the change in relative frequency is proportional to the product of fitness and frequency of a given type divided by the mean fitness.


class GRN:

    def __init__(self, Name, mRNA, Protein, Configuration,
                 Proportion, MutationRate, TranscriptionRate,
                 TranslationRate, DegradationRatemRNA,
                 DegradationRateProtein, DilutionRate,
                 TranscriptionThreshold, LogicGates, Sigmoid_k,
                 Leakage, f0, sys_input_ChIP, sys_PPI='', cis_impact=1):
        '''Parameters Initiation'''
        self.Name = Name
        self.mRNA = mRNA
        self.Protein = Protein
        self.cis_impact = cis_impact
        # Configuration.shape = (2,m,n); (0,m,n) for Activation and (1,m,n) for Repression
        self.Configuration = Configuration
        self.Proportion = Proportion
        self.MutationRate = MutationRate
        self.TranscriptionRate = TranscriptionRate
        self.TranslationRate = TranslationRate
        self.DegradationRatemRNA = DegradationRatemRNA
        self.DegradationRateProtein = DegradationRateProtein
        self.DilutionRate = DilutionRate
        # For each GRN, its shape should be (7000,7000), Genei is an activator of Genej.
        self.TranscriptionThreshold = TranscriptionThreshold
        self.LogicGates = LogicGates  # For each gene, it has 2 parameters
        self.Sigmoid_k = Sigmoid_k
        self.f0 = f0
        self.Leakage = Leakage
        self.Memo_mRNA = []
        self.sol = []

        if ((len(self.mRNA) + len(self.Protein)
            + len(self.TranscriptionRate) + len(self.TranslationRate)
            + len(self.DegradationRatemRNA) + len(self.DegradationRateProtein))
                / (6*len(self.mRNA)) != 1):

            raise Exception('len(mRNA) != len(Protein)')
        else:
            pass

        self.sys_input_ChIP = sys_input_ChIP
        self.PPI = sys_PPI

        return

    def ResetMemo_mRNA(self):
        self.Memo_mRNA = []
        self.sol = []
        return

    def UpdateProportion(self, NewProportion):
        '''Replace the old proportion by the new one'''
        self.Proportion = NewProportion
        if self.Proportion < 0 or self.Proportion > 1:
            raise Exception('Proportion Error!')
        return

    def GetDistance(self, WT_mRNA, transcriptionalprofilemax):
        '''Reads are equal'''
        # return(np.linalg.norm(self.mRNA - WT_mRNA))
        '''Genes are equal'''
        outdistance = 0.0
        for y in range(0, len(WT_mRNA)):
            outdistance = outdistance + \
                abs((self.mRNA[y]-WT_mRNA[y])/transcriptionalprofilemax[y])
        return outdistance

    def GetActivatorsRepressors(self):
        '''Return the indices of activators[0]/Repressors[1] for all genes'''
        arlist = [[] for i in range(0, self.Configuration.shape[0])]
        for i in range(0, self.Configuration.shape[0]):
            temp_list = []
            for j in range(0, len(self.f0)):
                temp_list.append(
                    np.where(self.Configuration[i][:, j] == 1)[0].tolist())
            arlist[i] = temp_list
        return arlist

    def Delta_mRNA(self, knockoutlist_in, TRM):
        '''Calculating the Delta mRNA'''
        ARList = self.GetActivatorsRepressors()
        scipystring = 'def update_mRNA_protein(t, y, knockoutlist):' + \
            '\n    delta_mRNA = [0]*int(len(y)/2)' + '\n    delta_protein = [0]*int(len(y)/2)' + \
            '\n    outlist = [0]*len(y)'  # Scipy
        for i in range(0, len(self.f0)):
            '''loop the mRNA list (For each gene)'''
            # Current Activators are ARList[0][i]; Current Repressors are ARList[1][i]

            '''Form the differential equations'''
            TempString_CA_sci = ''  # Scipy
            TempString_CR_sci = ''  # Scipy
            TempString_sci = ''    # Scipy

            '''Determine CA'''
            CA_complexes = {}
            for j in range(0, len(ARList[0][i])):
                if self.LogicGates[i][ARList[0][i][j]] in CA_complexes:
                    CA_complexes[self.LogicGates[i][ARList[0][i][j]]].append(ARList[0][i][j])
                else:
                    CA_complexes[self.LogicGates[i][ARList[0][i][j]]] = [ARList[0][i][j]]

            TempString_CSA_List = []
            for each_complex in CA_complexes:
                TempString_CA_complexes = ''
                for each_subunit in CA_complexes[each_complex]:
                    TempString_CA_complexes = 'ActivatorSigmoid({}, {}, {})*'.format(
                        'y[{}]'.format(len(self.f0)+each_subunit),
                        self.TranscriptionThreshold[each_subunit][i],
                        self.Sigmoid_k[each_subunit]) + TempString_CA_complexes  # Scipy
                TempString_CA_complexes = TempString_CA_complexes[:-1]
                TempString_CSA_List.append(TempString_CA_complexes)

            for each in TempString_CSA_List:
                TempString_CA_sci = '(1-{})*'.format(each)
            TempString_CA_sci = TempString_CA_sci[:-1]
            if TempString_CA_sci == '':
                TempString_CA_sci = '1'
            else:
                pass
            TempString_CA_sci = '(1-' + TempString_CA_sci + ')'

            '''Determine CR'''
            CR_complexes = {}
            for j in range(0, len(ARList[1][i])):
                if self.LogicGates[i][ARList[1][i][j]] in CR_complexes:
                    CR_complexes[self.LogicGates[i][ARList[1][i][j]]].append(ARList[1][i][j])
                else:
                    CR_complexes[self.LogicGates[i][ARList[1][i][j]]] = [ARList[1][i][j]]

            TempString_CSR_List = []
            for each_complex in CR_complexes:
                TempString_CR_complexes = ''
                for each_subunit in CR_complexes[each_complex]:
                    TempString_CR_complexes = 'ActivatorSigmoid({}, {}, {})*'.format(
                        'y[{}]'.format(len(self.f0)+each_subunit),
                        self.TranscriptionThreshold[each_subunit][i],
                        self.Sigmoid_k[each_subunit]) + TempString_CR_complexes  # Scipy
                TempString_CR_complexes = TempString_CR_complexes[:-1]
                TempString_CR_complexes = '(1-' + TempString_CR_complexes + ')'
                TempString_CSR_List.append(TempString_CR_complexes)
            for each in TempString_CSR_List:
                TempString_CR_sci = '{}*'.format(each)
            TempString_CR_sci = TempString_CR_sci[:-1]
            if TempString_CR_sci == '':
                TempString_CR_sci = '1'
            else:
                pass
            TempString_CR_sci = '(' + TempString_CR_sci + ')'

            TempString_sci = ('{}+({}-{})*('.format(self.Leakage[i],
                                                    self.TranscriptionRate[i],
                                                    self.Leakage[i])
                              + '{}+{}*(-(1-{})*(1-{}))+(1-{})*{}*{}'.format(
                self.f0[i],
                self.f0[i],
                TempString_CA_sci,
                TempString_CR_sci,
                self.f0[i],
                TempString_CA_sci,
                TempString_CR_sci)
                + ')-{}*y[{}]'.format(self.DegradationRatemRNA[i], i))  # Scipy

            proteinstring = '{}*y[{}] - ({}+{})*y[{}]'.format(self.TranslationRate[i],
                                                              i,
                                                              self.DegradationRateProtein[i],
                                                              self.DilutionRate[i],
                                                              i+len(self.f0))

            scipystring = (scipystring
                           + '\n    delta_mRNA[{}]='.format(i)
                           + TempString_sci
                           + '\n    delta_protein[{}]='.format(i)
                           + proteinstring)

        for i in range(0, len(self.f0)):
            scipystring = scipystring + ('\n    if (y[{}]<=0 and delta_mRNA[{}]<=0) or ({} in knockoutlist):'
                                         '\n        outlist[{}] = -y[{}]'
                                         '\n    elif {} in knockoutlist:'
                                         '\n        outlist[{}] = {}-y[{}]'
                                         '\n    else:'
                                         '\n        outlist[{}] = delta_mRNA[{}]'
                                         '\n    if y[{}]<=0 and delta_protein[{}]<=0:'
                                         '\n        outlist[{}]=-y[{}]'
                                         '\n    else:'
                                         '\n        outlist[{}]=delta_protein[{}]').format(
                i, i, i, i, i,
                -2-i, i, TRM[i], i, i,
                i, i+len(self.f0), i, i+len(self.f0), i+len(self.f0),
                i+len(self.f0), i)

        scipystring = scipystring + '\n    return outlist'
        return(scipystring)

    def Updatef0_1(self):
        '''Given a GRN and an attractor, return the f0 based on the steady-state assumption for each gene.'''
        ARList = self.GetActivatorsRepressors()
        f0_list = []
        for i in range(0, len(self.f0)):
            '''Form the differential equations'''
            TempString_CA_sci = ''  # Scipy
            TempString_CR_sci = ''  # Scipy

            '''Determine CA'''
            CA_complexes = {}
            for j in range(0, len(ARList[0][i])):
                if self.LogicGates[i][ARList[0][i][j]] in CA_complexes:
                    CA_complexes[self.LogicGates[i][ARList[0][i][j]]].append(ARList[0][i][j])
                else:
                    CA_complexes[self.LogicGates[i][ARList[0][i][j]]] = [ARList[0][i][j]]

            TempString_CSA_List = []
            for each_complex in CA_complexes:
                TempString_CA_complexes = ''
                for each_subunit in CA_complexes[each_complex]:
                    TempString_CA_complexes = 'ActivatorSigmoid({}, {}, {})*'.format(
                        self.Protein[each_subunit],
                        self.TranscriptionThreshold[each_subunit][i],
                        self.Sigmoid_k[each_subunit]) + TempString_CA_complexes  # Scipy
                TempString_CA_complexes = TempString_CA_complexes[:-1]
                TempString_CSA_List.append(TempString_CA_complexes)

            for each in TempString_CSA_List:
                TempString_CA_sci = '(1-{})*'.format(each)
            TempString_CA_sci = TempString_CA_sci[:-1]
            if TempString_CA_sci == '':
                TempString_CA_sci = '1'
            else:
                pass
            TempString_CA_sci = '(1-' + TempString_CA_sci + ')'

            '''Determine CR'''
            CR_complexes = {}
            for j in range(0, len(ARList[1][i])):
                if self.LogicGates[i][ARList[1][i][j]] in CR_complexes:
                    CR_complexes[self.LogicGates[i][ARList[1][i][j]]].append(ARList[1][i][j])
                else:
                    CR_complexes[self.LogicGates[i][ARList[1][i][j]]] = [ARList[1][i][j]]

            TempString_CSR_List = []
            for each_complex in CR_complexes:
                TempString_CR_complexes = ''
                for each_subunit in CR_complexes[each_complex]:
                    TempString_CR_complexes = 'ActivatorSigmoid({}, {}, {})*'.format(
                        self.Protein[each_subunit],
                        self.TranscriptionThreshold[each_subunit][i],
                        self.Sigmoid_k[each_subunit]) + TempString_CR_complexes  # Scipy
                TempString_CR_complexes = TempString_CR_complexes[:-1]
                TempString_CR_complexes = '(1-' + TempString_CR_complexes + ')'
                TempString_CSR_List.append(TempString_CR_complexes)

            for each in TempString_CSR_List:
                TempString_CR_sci = '{}*'.format(each)
            TempString_CR_sci = TempString_CR_sci[:-1]
            if TempString_CR_sci == '':
                TempString_CR_sci = '1'
            else:
                pass
            TempString_CR_sci = '(' + TempString_CR_sci + ')'
            CA = eval(TempString_CA_sci)
            CR = eval(TempString_CR_sci)
            temp_value = np.round((((self.DegradationRatemRNA[i]*self.mRNA[i]-self.Leakage[i])/(
                self.TranscriptionRate[i]-self.Leakage[i]))-CA*CR)/(CA+CR-2*CA*CR), 5)
            if abs(CA+CR-2*CA*CR) < 0.01:
                f0_list.append('mask')
            elif math.isnan(temp_value):
                f0_list.append(np.nan)
            elif temp_value < 0:
                f0_list.append(np.nan)
            elif temp_value > 1:
                f0_list.append(np.nan)
            else:
                f0_list.append(temp_value)
        return(f0_list)

    def Updatef0_2(self, newf0):
        self.f0 = newf0
        return

    def GenerateMutation(self, globalmutationdirection, samplesize, globalmutationdirection_LG):
        '''Mutations occur upon 1.Configuration; 2.LogicGates'''

        '''1.Mutation upon Configuration'''
        self.Configuration = HammingMutation(
            self.sys_input_ChIP, self.MutationRate, self.Configuration,
            globalmutationdirection, samplesize, 0.99)

        '''2.Mutation upon LogicGates'''
        self.LogicGates = Mutation_LG(self.MutationRate, self.LogicGates, self.GetActivatorsRepressors(), self.PPI)

        return

    def SetmRNA(self, NewmRNAList, knockout_list, Overexpression):
        #NewmRNAList = NewmRNAList.tolist()
        self.mRNA = copy.deepcopy(NewmRNAList)
        for i in range(0, len(knockout_list)):
            if knockout_list[i] >= 0:
                self.mRNA[knockout_list[i]] = 0
            else:
                self.mRNA[-knockout_list[i]-2] = Overexpression[i]
        return

    def SetmRNA_AntiSelfAct(self, NewmRNAList, knockout_list, TMax, overexpression, WTTP_, sysi):
        '''Update mRNA in a way that accounts self activation'''
        #NewmRNAList = NewmRNAList.tolist()
        self.mRNA = [copy.deepcopy(NewmRNAList)]
        temp_mRNA = copy.deepcopy(NewmRNAList)
        to_append = 0
        for i in range(0, len(self.mRNA[0])):

            if ((self.Configuration[0][i][i] == 1 and self.Configuration[1][i][i] == 0) and (self.LogicGates[i][0] == 0 or np.count_nonzero(self.Configuration[0][:, i]) < 2)):
                Comparison_List = []
                Comp_mRNA = copy.deepcopy(WTTP_[str(sysi)][1])
                Comp_mRNA[i] = max(TMax[i], overexpression[i])
                for j in range(0, len(WTTP_)):
                    Comparison_List.append(GetAttractorDistance(
                        Comp_mRNA, WTTP_[str(j)][1], TMax)/len(Comp_mRNA))
                if min(Comparison_List) < 0.1:
                    pass
                else:
                    temp_mRNA[i] = TMax[i]
                    to_append = 1
                    # self.mRNA.append(copy.deepcopy(temp_mRNA))
            else:
                pass

        if to_append == 1:
            self.mRNA.append(copy.deepcopy(temp_mRNA))
        else:
            pass
        for z in range(0, len(self.mRNA)):
            for i in range(0, len(knockout_list)):
                if knockout_list[i] >= 0:
                    self.mRNA[z][knockout_list[i]] = 0
                else:
                    self.mRNA[z][-knockout_list[i]-2] = TMax[i]
        return

    def SetMemo(self, NewMemo):
        self.Memo_mRNA = NewMemo
        return

    def SetProtein(self, NewProteinList, knockout_list, Overexpression):
        #NewProteinList = NewProteinList.tolist()
        self.Protein = copy.deepcopy(NewProteinList)
        for i in range(0, len(knockout_list)):
            if knockout_list[i] >= 0:
                self.Protein[knockout_list[i]] = 0
            else:
                j = int(-knockout_list[i]-2)
                self.Protein[j] = (self.TranslationRate[j]*Overexpression[j]) / \
                    (self.DegradationRateProtein[j]+self.DilutionRate[j])
        return

    def SetProtein_AntiSelfAct(self, knockout_list, Overexpression):
        '''Update protein values in a way that accounts for self-activation'''
        self.Protein = []
        for i in range(0, len(self.mRNA)):
            NewProtein_ = np.array(
                [0 for x in range(0, len(self.mRNA[0]))], dtype=float)
            for j in range(0, len(NewProtein_)):
                NewProtein_[j] = (self.TranslationRate[j]*self.mRNA[i][j]) / \
                    (self.DegradationRateProtein[j]+self.DilutionRate[j])
            self.Protein.append(NewProtein_)

        for j in range(0, len(self.Protein)):
            for i in range(0, len(knockout_list)):
                if knockout_list[i] >= 0:
                    self.Protein[j][knockout_list[i]] = 0
                else:
                    z = int(-knockout_list[i]-2)
                    self.Protein[j][z] = (self.TranslationRate[z]*Overexpression[z])/(
                        self.DegradationRateProtein[z]+self.DilutionRate[z])
        return

    def SetMutationRate(self, NewMutationRate):
        self.MutationRate = copy.deepcopy(NewMutationRate)
        return

    def SetTranslationRate(self, NewTranslationRate):
        self.TranslationRate = copy.deepcopy(NewTranslationRate)
        return

    def Update_Transcription_Rate(self, TranscriptionRate):
        self.TranscriptionRate = TranscriptionRate
        return

    def Update_Promoter_Strength(self, which_i, sys_promoter_strengths, sys_mRNA_elongation_rate, sys_gene_length):
        TranscriptionRate = sys_promoter_strengths[which_i]*60*sys_mRNA_elongation_rate/sys_gene_length.tolist()
        self.TranscriptionRate = TranscriptionRate
        return

