# -*- coding: utf-8 -*-
"""
# Copyright 2013,2014,2015 Georgia Institute of Technology.
#
# Licensed under the Apache License, Version 2.0 (the “License”); you may not use
# this file except in compliance with the License. You may obtain a copy of the
# License at http://www.apache.org/licenses/LICENSE-2.0. Unless required by
# applicable law or agreed to in writing, software distributed under the License
# is distributed on as “AS IS” BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied.  See the License for the specific language
# governing permissions and limitations under the License.

Created on Tue Jan 13 13:24:07 2015

The function computes varies kinds of information-based clustering quality measures.
The function currently supports:

    MI:     mutual information
    NMI:    normalized mutual information
    ID:     information distance
    NID:    normalized information distance
    VI:     variation of information
    NVI:    normalized variation of information
"""

import numpy as np
import sys
import getopt

# define constants
macheps = np.finfo(float).eps

class contingency_table:
    def __init__(self, clusterLabels1, clusterLabels2):
        self.clusterLabels1 = clusterLabels1
        self.clusterLabels2 = clusterLabels2
    # end ctor
        
    def contingency(self): 
        # create a contingency table    
        clusterLabels1 = self.clusterLabels1
        clusterLabels2 = self.clusterLabels2
        
        mem1ShapeMin = np.min(np.shape(clusterLabels1))
        mem2ShapeMin = np.min(np.shape(clusterLabels2))
        
        if (mem1ShapeMin > 1) or (mem2ShapeMin > 1):
            sys.exit('Contingency: Requires two vectors (nx1 or 1xn) arguments')

        # add an extra clumn and an extra row for cluster label being 0
        contingencyTable = np.zeros(shape=(np.max(clusterLabels1)+1,np.max(clusterLabels2)+1))
    
        upperLoopIdx = np.arange(np.size(clusterLabels1))
        
        for i in upperLoopIdx: 
            mem1Value = clusterLabels1[0,i]
            mem2Value = clusterLabels2[0,i] 
            contingencyTable[mem1Value, mem2Value] = contingencyTable[mem1Value, mem2Value] + 1
            
        # determine row(s) and/or column(s) of the contingency table that do 
        # not contain any samples which could be removed
        rowsToKeep, colsToKeep = self.missing_columns_rows()
        
        # remove those row(s) and/or column(s) from the contingency table
        contingencyTableOut = self.remove_missing_columns_rows(contingencyTable, \
            rowsToKeep, colsToKeep)

        return contingencyTableOut
    # end contingency() 
    
    ####### Helper methods #######
    def ismember(self, clusterRange, clusterLabels):
        # provides an ismember method as in Matlab(tm)
        bind = {}
    
        for i, elt in enumerate(clusterLabels):
            if elt not in bind:
                bind[elt] = i
            
        memberIndices = [bind.get(itm, -1) for itm in clusterRange] 
    
        return memberIndices
    # end ismember()
        
    def missing_columns_rows(self):
        # determine row(s) and/or columns(s) of the contingency table that are
        # to be kept    
        clusterLabels1 = list(np.squeeze(np.asarray(self.clusterLabels1)))
        clusterLabels2 = list(np.squeeze(np.asarray(self.clusterLabels2)))

        maxClusterLabel1 = int(np.max(clusterLabels1))
        clusterLabel1Range = range(maxClusterLabel1+1)
    
        rowMask = self.ismember(clusterLabel1Range, clusterLabels1)
    
        # add 1 to list; turns -1s into 0
        rowMaskMod = [x+1 for x in rowMask]
        # convert list to boolean array. now 0s are False; other vals True
        rowsToKeep = np.asarray(rowMaskMod, dtype=bool)
    
        maxClusterLabel2 = int(np.max(clusterLabels2))
        clusterLabel2Range = range(maxClusterLabel2+1)
        colMask = self.ismember(clusterLabel2Range, clusterLabels2)
    
        # add 1 to list; turns -1s into 0
        colMaskMod = [x+1 for x in colMask]
        # convert list to boolean array. now 0s are False; other vals True
        colsToKeep = np.asarray(colMaskMod, dtype=bool)

        return rowsToKeep, colsToKeep
    # end missing_columns_rows()
        
    def remove_missing_columns_rows(self, contingencyTable, rowsToKeep, colsToKeep):
        # remove rows(s) and/or column(s) of the contingency table that do not
        # contain any numbers. corresponds to missing labels
    
        contingencyTable1 = contingencyTable[rowsToKeep,:]
        contingencyTable1 = contingencyTable1[:,colsToKeep]
    
        return contingencyTable1
    # end remove_missing_columns_rows()
        
# end class contingency_table

class entropy:
    # constructor for the entropy class
    # 2015/02/22: all entropies have class members
    def __init__(self, clusterLabels1, clusterLabels2):
        # compute the contingency table based on clusterLabels1 and clusterLabels2
        self.clusterLabels1 = clusterLabels1
        self.clusterLabels2 = clusterLabels2
        cContingency = contingency_table(clusterLabels1, clusterLabels2)
        contingencyTable = cContingency.contingency()
        
        # compute row and column sums, and determine total number of samples 
        # from contingencyTable
        numSamples = np.sum(np.sum(contingencyTable))
        numRows, numCols = np.shape(contingencyTable)
    
        # compute the row sums and col sums of contingency table
        if numCols > 1:
            rowSums = np.sum(contingencyTable.T, axis=0)
        else:
            rowSums = contingencyTable[:,0]
        
        if numRows > 1:
            colSums = np.sum(contingencyTable, axis=0)
        else:
            colSums = contingencyTable[0,:]
            
        self.contingencyTable = contingencyTable
        self.numSamples = numSamples
        self.numRows = numRows
        self.numCols = numCols
        self.rowSums = rowSums
        self.colSums = colSums

        # for accessor to all entropies
        self.entropy = 0.0
        self.entropy1 = 0.0
        self.entropy2 = 0.0
        self.jointEntropy = 0.0
        self.conditionalEntropy12 = 0.0
        self.conditionalEntropy21 = 0.0
    # end ctor
        
    def compute_entropy(self, clusterResult=1): 
        # for a cluster result U compute the entropy H(U)    
        # compute probability mass function and its log
        if int(clusterResult) == 1:
            # compute probability mass function of the first clustering result
            probMass = self.rowSums / self.numSamples
        else:
            # compute probability mass function of the second clustering result
            probMass = self.colSums / self.numSamples
        
        logProbMass = np.log(probMass) 
        
        # compute entropy
        self.entropy = np.dot(-probMass, logProbMass.T) 
    
        return self.entropy
    # end compute_entropy()
        
    def compute_entropies(self):
        # for clusterings results U and V compute the entropy H(U) and H(V)
        # computedEntropy1 for the first clustering result U, and computedEntropy2
        # the second
    
        # compute probability mass functions and their logs
        probMass1 = self.rowSums / self.numSamples
        logProbMass1 = np.log(probMass1) 
        probMass2 = self.colSums / self.numSamples
        logProbMass2 = np.log(probMass2) 

        # compute entropies
        self.entropy1 = np.dot(-probMass1, logProbMass1.T) 
        self.entropy2 = np.dot(-probMass2, logProbMass2.T)
    
        return self.entropy1, self.entropy2
    # end compute_entropies()
        
    def joint_entropy(self):
        # For two clusterings U and V compute H(U,V) 
    
        jointEntropy = 0.0
        for i in np.arange(self.numRows):
            for j in np.arange(self.numCols):
                if self.contingencyTable[i,j] > 0:
                    jointEntropy = jointEntropy + -1.0*(self.contingencyTable[i,j] \
                        * np.log(self.contingencyTable[i,j] / self.numSamples))

        self.jointEntropy = jointEntropy / self.numSamples
        
        return self.jointEntropy
    # end joint_entropy()
    
    def conditional_entropies(self):
        # For two clusterings U and V, compute H(U|V) (entropy of U conditioned on V)
        # U being the first clustering result, and V the second, if clusterNumToCondition = 1,
        # else H(V|U) would be computed. Note that H(U|V) != H(V|U), necessarily, i.e.
        # not symmetric.
        conditionalEntropy12 = 0.0
        for i in np.arange(self.numRows):
            for j in np.arange(self.numCols):
                if self.contingencyTable[i,j] > 0:
                    conditionalEntropy12 = conditionalEntropy12 + -1.0*(self.contingencyTable[i,j] \
                        * np.log(self.contingencyTable[i,j] / self.colSums[j]))

        self.conditionalEntropy12 = conditionalEntropy12 / self.numSamples

        conditionalEntropy21 = 0.0
        # contingency table has to be transposed 
        contingencyTable = self.contingencyTable.T
        numRows = self.numCols
        numCols = self.numRows
        for i in np.arange(numRows):
            for j in np.arange(numCols):
                if contingencyTable[i,j] > 0:
                    conditionalEntropy21 = conditionalEntropy21 + -1.0*(contingencyTable[i,j] \
                        * np.log(contingencyTable[i,j] / self.rowSums[j]))

        self.conditionalEntropy21 = conditionalEntropy21 / self.numSamples
        
        return self.conditionalEntropy12, self.conditionalEntropy21
    # end conditional_entropy()

    def get_entropies(self):
        # return all the entropies    
        allEntropies = {}
    
        allEntropies =  {'entropy':self.entropy,
                         'entropy1': self.entropy1,
                         'entropy2': self.entropy2,
                         'joint': self.jointEntropy,
                         'conditional1': self.conditionalEntropy12,
                         'conditional2': self.conditionalEntropy21} 
        return allEntropies
    # end get_entropies()        
# end class entropy

class measures:
    # constructor for the entropy class
    # 2015/02/22: added method to return "pointer" to entropy instance
    def __init__(self, clusterLabels1, clusterLabels2):
        # Inputs:
        # clusterLabels1: integer values assigned to each cluster of the first clustering result
        # clusterLabels2: integer values assigned to each cluster of the second clustering result 
        self.clusterLabels1 = np.asmatrix(clusterLabels1) 
        self.clusterLabels2 = np.asmatrix(clusterLabels2) 
        
        # call the class entropy to initialize necessary entropies
        self.cEntropy = entropy(self.clusterLabels1, self.clusterLabels2)
        entropy1, entropy2 = self.cEntropy.compute_entropies()
        jointEntropy = self.cEntropy.joint_entropy()
        ce1, ce2 = self.cEntropy.conditional_entropies()
        
        self.entropy1 = entropy1
        self.entropy2 = entropy2
        self.jointEntropy = jointEntropy
        self.conditionalEntropy1 = ce1
        self.conditionalEntropy2 = ce2
    # end ctor
            
    def getEntropyPtr(self):
        return self.cEntropy
            
    def getAllMeasures(self):
        dAllMeas =  {'MI': self.mutual_information(),
                         'NMI': self.normalized_mutual_information(),
                         'ID': self.information_distance(),
                         'NID': self.normalized_information_distance(),
                         'VI': self.variation_information(),
                         'NVI': self.normalized_variation_information()} 
    
        #print "************* exiting method test_measures ****************\n"
        return dAllMeas

    def mutual_information(self):
        # compute mutual information (MI) measure.
        # A higher value of MI indicates that the two sets of cluster labels  
        # are nonrandomly associated to each other. The higher the MI, the more 
        # useful the information in one cluster labels set helps to predict the
        # cluster labels in another and vice-versa; its range is 
        # [0, min(entropy1, entropy2)]. Note that only one
        # of the twe computed conditional entropies are needed to comute MI
        # since the 'intersection' of information remains the same.    
        mutualInformation = self.entropy1 - self.conditionalEntropy1
        mutualInformation = np.fabs(mutualInformation) # fix display of zero for special case         
        return mutualInformation
    # end mutual_information()
        
    def normalized_mutual_information(self):
        # compute normalized mutual information (MI) measure, its range is [0,1].
        # As it is closer to 1, the two cluster labels set are more nonrandomly
        # associated with each other
        
        # compute mutual information
        mutualInformation = self.mutual_information()
        minEntropy = np.min([self.entropy1, self.entropy2])
        
        # test for near machine epsilon value; use alternative computation
        if  (np.fabs(minEntropy) + 1.0) > 1.0:
            # we use the square root form
            normalizedMutualInfo = mutualInformation / np.sqrt(self.entropy1*self.entropy2)
            normalizedMutualInfo = np.fabs(normalizedMutualInfo) # fix display of zero for special case         
        elif ((np.fabs(self.entropy1) + 1.0) > 1.0) or ((np.fabs(self.entropy2) + 1.0) > 1.0):
            warningString = """***Warning***: One of the entropies computed is close to 0,
                normalized mutual information computed would be divided by a
                very small number instead of the entropy. Alternative
                computation being used for MI.\n"""
            print(warningString)
                
            # Use an alternative formulation of NMI
            normalizedMutualInfo = mutualInformation / (self.entropy1 + self.entropy2)
            normalizedMutualInfo = np.fabs(normalizedMutualInfo) # fix display of zero for special case         
        else:
            warningString = """***Warning***: The entropies computed are close to 0. Please check 
            the inputs; setting the result to undefined\n\n"""
            print(warningString)
            normalizedMutualInfo = 'undefined'

        return normalizedMutualInfo
    # end normalized_mutual_information()
        
    def information_distance(self):
        # compute information distance measure, its range is [0, log(number of samples)].
        # A value closer to 0 means the distance in terms of information between
        # the two sets of clusters are closer or more similar
    
        # compute mutual information
        mutualInformation = self.mutual_information()
        
        informationDistance = np.max([self.entropy1, self.entropy2]) - mutualInformation
        
        informationDistance = np.fabs(informationDistance) # fix display of zero for special case         
        return informationDistance
    # end information_distance()
    
    def normalized_information_distance(self):
        # compute normalized information distance measure, its range is [0,1].
        # A value closer to 0 means the NID is lower between the 
        # two sets of cluster labels; the two clustering are more similar
    
        # compute mutual information
        mutualInformation = self.mutual_information()        
        maxEntropy = np.max([self.entropy1, self.entropy2])
        
        if (np.fabs(maxEntropy) + 1.0) > 1.0:
            normalizedInfoDistance = 1.0 - mutualInformation / maxEntropy
        else:
            warningString = """***Warning***: The entropies computed are close to 0. Please check 
            the inputs; setting the result to undefined\n"""
            print(warningString)
            normalizedInfoDistance = 'undefined'
            
        return normalizedInfoDistance
    # end normalized_information_distance()
        
    def variation_information(self):
        # computes the variation of information measure; the value is 0 iff two sets of
        # clustering labels are the same, its range is [0, log(number of samples)].
        # A value closer to 0 means the VI is lower between the 
        # two sets of cluster labels; the two clustering are more similar
    
        # compute mutual information
        mutualInformation = self.mutual_information()
        variationInformation = self.jointEntropy - mutualInformation
            
        return variationInformation
    # end variation_information()
    
    def normalized_variation_information(self):
        # compute normalized variation of information measure, its range is [0, 1].
        # A value closer to 0 means the NVI is lower between the 
        # two sets of cluster labels; the two clustering are more similar
    
        # compute mutual information
        mutualInformation = self.mutual_information()
        
        if (np.fabs(self.jointEntropy) + 1.0) > 1.0:
            normalizedVariationInfo = 1.0 - (mutualInformation / self.jointEntropy)
        else:
            warningString = """***Warning***: The entropies computed are close to 0. Please check 
            the inputs; setting the result to undefined\n\n"""
            print(warningString)
            normalizedVariationInfo = 'undefined'
                
            # No alternative computation that approximates this distance metric.
         
        return normalizedVariationInfo
    # end normalized_variation_information()            
# end class measures
        
class outlier_adjustments:
    def __init__(self, clusterLabels1, clusterLabels2, outlierClusterLabel=-1):
        
        self.clusterLabels1 = clusterLabels1
        self.clusterLabels2 = clusterLabels2
        self.outlierClusterLabel = outlierClusterLabel
        
    def modify_outliers_cluster_label(self, clusterLabelsToAdjust=1):
        # the function change a set of cluster labels for the outliers from 
        # outlierClusterLabel to the maximum cluster label number + 1 so the 
        # contingency table created later would account those outliers in the 
        # last row and/or column of the contingency table. The outliers would 
        # most likely affect the quality measures
    
        if int(clusterLabelsToAdjust) == 1:
            # cluster labels to adjust are the first set of cluster labels
            clusterLabels = self.clusterLabels1
        else:
            # cluster labels to adjust are the second set of cluster labels
            clusterLabels = self.clusterLabels2

        maxClusterLabel = np.max(clusterLabels)
    
        # find the document indices which have clusterLabelsToAdjust cluster 
        # labels (outliers)
        outlierIndices = np.where(clusterLabels == self.outlierClusterLabel)[0]
    
        # change outliers cluster label from clusterLabelsToAdjust to a large 
        # number incremented by 1 from maxClusterLabel
        clusterLabels[outlierIndices] = maxClusterLabel + 1
    
        return clusterLabels
    # end modify_outliers_cluster_label()
    
    def remove_outliers(self):
        # the function remove documents that are considered as outliers from  
        # both sets of cluster labels. If a document is considered as an  
        # outlier in clusterLabels1 but not in clusterLabels2, the document 
        # would still be removed in clusterLabels2, and vice-versa
    
        clusterLabels1 = self.clusterLabels1
        clusterLabels2 = self.clusterLabels2

        # find documents indices that are not outliers
        keepIndices1 = np.where(clusterLabels1 != self.outlierClusterLabel)[0]
        keepIndices2 = np.where(clusterLabels2 != self.outlierClusterLabel)[0]
    
        # find document indices that are not outliers in both of the cluster sets
        keepIndices = set(keepIndices1) & set(keepIndices2)
    
        # turn those keepIndices by to an array
        keepIndices = np.array(list(keepIndices), dtype=int)
    
        # keep those documents with indices in keepIndices only
        clusterLabels1Out = clusterLabels1[keepIndices]
        clusterLabels2Out = clusterLabels2[keepIndices]
    
        return clusterLabels1Out, clusterLabels2Out
    # end remove_outliers()
# end class outlier_adjustments

def load_assignment_files(filesDir, fileName1, fileName2):
    # the function reads in two files to obtain two sets of cluster assignment
    # labels    
    import csv
        
    # extract first row in the first file for the first set of cluster labels
    stop = False
    #print 'fileName1: ', fileName1
    thePath = filesDir + '/' + fileName1
    print 'fileName1: ', thePath
    try:
        with open(thePath, 'rb') as f1:
            reader = csv.reader(f1)
            for row in reader:
                stop = True
                clusterLabels1 = row
                if stop == True:
                    break
    except:
        sys.exit('Error: Could not find the file specified, please check')
                    
    #print 'fileName2: ', fileName2
    thePath = filesDir + '/' + fileName2
    print 'fileName2: ', thePath
    # extract first row in the second file for the second set of cluster labels        
    stop = False
    try:
        with open(thePath, 'rb') as f2:
            reader = csv.reader(f2)
            for row in reader:
                stop = True
                clusterLabels2 = row
                if stop == True:
                    break
    except:
        sys.exit('Error: Could not find the file specified, please check')
    
    # ensure the cluster labels could be converted to integers
    try:
        clusterLabels1 = np.array(clusterLabels1, dtype=np.int)
        clusterLabels2 = np.array(clusterLabels2, dtype=np.int)
    except:
        sys.exit('Error: The cluster labels in the file(s) may contain non-integer values, please check')   
    
    # ensure the two sets of cluster labels have the same length
    if len(clusterLabels1) != len(clusterLabels2):
        sys.exit('Error: The lengths of the two sets of cluster labels differ, please check')
    
    return clusterLabels1, clusterLabels2
# end load_assignment_files

def load_assignment_dir(filesDir):
    # the function reads in all files in a directory to obtain a list of cluster assignment
    # labels. the files have names like 'news20_assignments_32_run1.csv'
    
    # local imports
    import os
        
    # extract first row in the first file for the first set of cluster labels
    lFiles = []
    
    for root, dirs, files in os.walk(filesDir):
        for file in files:
            if file.endswith('.csv'):
                lFiles.append(file)

    return lFiles
# end  load_assignment_dir()
    
def test_measures(toLoadFiles, filesDir, fileName1, fileName2, clusterLabels1, clusterLabels2):
    # test function to test quality_measures.py; bypasses command line if desired
    # run with:
#        clusterLabels1 = [0,1,3,0,6,5,4,14,10,23,22,11,17,20]
#        clusterLabels2 = [1,0,5,4,3,0,9,19,23,26,21,19,20,28]

    
    if toLoadFiles == True:
        # load two sets of cluster labels for testing
        clusterLabels1, clusterLabels2 = load_assignment_files(filesDir, fileName1, fileName2)
        
        # let's test the quality measures by removing the outliers
        cOutlierAdjust = outlier_adjustments(clusterLabels1, clusterLabels2, \
            outlierClusterLabel=-1)
        clusterLabels11, clusterLabels22 = cOutlierAdjust.remove_outliers()

        # now test the quality measures by keeping the outliers, but adjust the
        # outlier cluster label
        clusterLabels1 = cOutlierAdjust.modify_outliers_cluster_label(clusterLabelsToAdjust=1)
        clusterLabels2 = cOutlierAdjust.modify_outliers_cluster_label(clusterLabelsToAdjust=2)
        
        # compute mutual information and all other supporting measures with outliers removed
        cMeasures1 = measures(clusterLabels11, clusterLabels22, clusterNumToCondition=1)
        mutualInformation = cMeasures1.mutual_information()
        informationDistance = cMeasures1.information_distance()
        variationInformation = cMeasures1.variation_information()
        normalizedMutualInfo = cMeasures1.normalized_mutual_information()
        normalizedInfoDistance = cMeasures1.normalized_information_distance()
        normalizedVariationInfo = cMeasures1.normalized_variation_information() 

        #print "************* exiting method test_measures ****************\n"
        dAllMeas =  {'MI': mutualInformation,
                         'NMI': normalizedMutualInfo,
                         'ID': informationDistance,
                         'NID': normalizedInfoDistance,
                         'VI': variationInformation,
                         'NVI': normalizedVariationInfo} 

        return dAllMeas
    # end if toLoadFiles
    
    # compute mutual information and all other supporting measures
    cMeasures = measures(clusterLabels1, clusterLabels2)
    mutualInformation = cMeasures.mutual_information()
    informationDistance = cMeasures.information_distance()
    variationInformation = cMeasures.variation_information()
    normalizedMutualInfo = cMeasures.normalized_mutual_information()
    normalizedInfoDistance = cMeasures.normalized_information_distance()
    normalizedVariationInfo = cMeasures.normalized_variation_information()

    dAllMeas =  {'MI': mutualInformation,
                     'NMI': normalizedMutualInfo,
                     'ID': informationDistance,
                     'NID': normalizedInfoDistance,
                     'VI': variationInformation,
                     'NVI': normalizedVariationInfo} 

    #print "************* exiting method test_measures ****************\n"
    return dAllMeas
    
# end test_measures()   
            
def quality_measures(metricId, bNoOutliers, filesDir, fileName1, fileName2):
    # the function calls to compute necessary clustering quality measures
    # 2015/02/22: Now returns all entropies

    # check values in the metricId first to ensure they are appropriate
    metricIdSet = set(metricId)
    supportingValueSet = set([1,2,3,4,5,6])
    if metricIdSet.issubset(supportingValueSet) == False:
        sys.exit('Error: The metricId input contains invalid value(s), please check')

    metricId = np.array(metricId, dtype=int)
    
    # load two sets of cluster labels
    clusterLabels1, clusterLabels2 = load_assignment_files(filesDir, fileName1, fileName2)
    
    # determine if there is variations in clusterLabels1 and clusterLabels2
    if np.var(clusterLabels1) < macheps or np.var(clusterLabels2) < macheps:
        warningString = """***Warning***: One or both of the cluster label sets contain only a
            single cluster, the quality measures may be inaccurate, please check\n\n"""
        print(warningString)
        
    # let's test the quality measures by removing the outliers
    cOutlierAdjust = outlier_adjustments(clusterLabels1, clusterLabels2, \
        outlierClusterLabel=-1)
            
    if bNoOutliers:
        clusterLabels1, clusterLabels2 = cOutlierAdjust.remove_outliers()
    else:
        print "outliers will be used for results.\n"
        clusterLabels1 = cOutlierAdjust.modify_outliers_cluster_label(1)
        clusterLabels2 = cOutlierAdjust.modify_outliers_cluster_label(2)
    
    if any(clusterLabels1 < np.int(0)) or any(clusterLabels2 < np.int(0)):
        sys.exit('Error: The cluster labels contain negative integers, please check')
    
    # gather all the entropies
    cMeasures = measures(clusterLabels1, clusterLabels2) # object to access measure methods
    entPtr = cMeasures.getEntropyPtr() # get pointer to entropy object
    allEntropies = entPtr.get_entropies() # accessor for all entropies
    
    # compute selected quality measures
    measuresAll = {}
    if any(metricId == 1):
        # compute mutual information measure
        mutualInformation = cMeasures.mutual_information()
        measuresAll['mutualInformation'] = mutualInformation
    if any(metricId == 2):
        # compute information distance measure
        informationDistance = cMeasures.information_distance()
        measuresAll['informationDistance'] = informationDistance
    if any(metricId == 3):
        # compute variation information measure
        variationInformation = cMeasures.variation_information()
        measuresAll['variationInformation'] = variationInformation
    if any(metricId == 4):
        # compute normalized mutual information measure
        normalizedMutualInfo = cMeasures.normalized_mutual_information()
        measuresAll['normalizedMutualInfo'] = normalizedMutualInfo
    if any(metricId == 5):
        # compute normalized information distance measure
        normalizedInfoDistance = cMeasures.normalized_information_distance()
        measuresAll['normalizedInfoDistance'] = normalizedInfoDistance
    if any(metricId == 6):
        # compute normalized variation information measure
        normalizedVariationInfo = cMeasures.normalized_variation_information()
        measuresAll['normalizedVariationInfo'] = normalizedVariationInfo
        
    return measuresAll, allEntropies
    
# end quality_measures()

###############################################################################
# Command line documentation utils
# Put performance measures and visualization code here
def usage_metrics():
    strUsage ="""
        usage: python quality_measures.py [-h or --help] ...
        
        Clustering quality measures code developed at the Georgia Institute of Technology.
        
        arguments:
            -d  directory path    REQUIRED: Use this directory to access files
            
            -f filename1          REQUIRED: Filename of the first set of cluster labels

            -g filename2          REQUIRED: Filename of the second set of cluster labels
            
            -m metric             Metric to use. MI,NMI(default),ID,VI,NMI,NID,NVI,all
            
            -r <boolean>          Remove outliers: False or True(default)
            
            -e <boolean>          Show entropies: False(default) or True
            
        optional arguments:
            -h, --help            Show this help message and exit
    """
    print strUsage
# end usage_framwork()
    
###############################################################################
# main()
# use to run tests or run application code                
def main(argv):
    # parse command line
    metric = None
    
    try:
        opts, args = getopt.getopt(argv, "hd:f:g:m:r:e:",['help','directory','filename1','filename2','metric, entropies'])
    except:
        usage_metrics()
        print 'argv is:', argv
        sys.exit(2)
    # end try
        
    # defalut values:
    bNoOutliers = True
    bShowEntropies = False
    metric = 'NMI'
    
    if not opts:
        print "No options provided. Using '-h' ", opts
        usage_metrics()
        return
        
    # check for REQUIRED options
    d_exists = False; f_exists = False; g_exists = False
    
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage_metrics()
            sys.exit()
        elif opt in ("-d", "--directory"):
            filesDir = arg
            d_exists = True # REQUIRED
        elif opt in ("-f", "--filename1"):
            filename1 = arg
            f_exists = True # REQUIRED
        elif opt in ("-g", "--filename2"):
            filename2 = arg
            g_exists = True # REQUIRED
        elif opt in ("-m", "--metric"):
            metric = arg
        elif opt in ("-r", "--remove"):
            #bNoOutliers = arg
            if arg == 'True': bNoOutliers = True
            if arg == 'False': bNoOutliers = False
        elif opt in ("-e", "--entropies"):
            #bNoOutliers = arg
            if arg == 'True': bShowEntropies = True
            if arg == 'False': bShowEntropies = False
        else:
            print "The command line option %s is not supported" % opt
            usage_metrics()
    # end for opts
        
    if not d_exists:
        print "****** A directory path is a REQUIRED option ******"
        usage_metrics()
        sys.exit()
        
    if not f_exists:
        print "****** A filename for the first set of cluster labels is a REQUIRED option ******"
        usage_metrics()
        sys.exit()
        
    if not g_exists:
        print "****** A filename for the second set of cluster labels is a REQUIRED option ******"
        usage_metrics()
        sys.exit()
        
    # end for opt
    # dMetricid = {1:'mutualInformation', 2:'informationDistance',3:'variationInformation',\
    #             4:'normalizedMutualInfo', 5:'normalizedInfoDistance',6:'normalizedVariationInfo'}
    metric = metric.upper() # user could lower of upper case; both should be valid
    metricId = []
    
    if not metric:
        metricId = [4]
        print "No metric spescified; computing NMI"
    elif metric == 'MI':
        metricId = [1]
    elif metric == 'ID':        
        metricId = [2]
    elif metric == 'VI':        
        metricId = [3]
    elif metric == 'NMI':        
        metricId = [4]
    elif metric == 'NID':        
        metricId = [5]
    elif metric == 'NVI':        
        metricId = [6]
    elif metric == 'all':        
        metricId = [1, 2, 3, 4, 5, 6]
    else:
        print "Metric %s not found" % metric

    print '\n'
    print "metric: ", metric
    #metricId = [1, 2, 3, 4, 5, 6]
        
    measuresAll, allEntropies = quality_measures(metricId, bNoOutliers, filesDir, filename1, filename2)
    print '\n'  
    if metric == 'MI': print 'Mutual Information: ', measuresAll['mutualInformation']
    if metric == 'NMI': print 'Normalized Mutual Information: ', measuresAll['normalizedMutualInfo']
    if metric == 'ID': print 'Information Distance: ', measuresAll['informationDistance']
    if metric == 'NID': print 'Normalized Information Distance: ', measuresAll['normalizedInfoDistance']
    if metric == 'VI': print 'Variation of Information: ', measuresAll['variationInformation']
    if metric == 'NVI': print 'Normalized Variation of Information: ', measuresAll['normalizedVariationInfo']
    if metric == 'all':
        print 'Mutual Information:                  ', measuresAll['mutualInformation']
        print 'Normalized Mutual Information:       ', measuresAll['normalizedMutualInfo']
        print 'Information Distance:                ', measuresAll['informationDistance']
        print 'Normalized Information Distance:     ', measuresAll['normalizedInfoDistance']
        print 'Variation of Information:            ', measuresAll['variationInformation']
        print 'Normalized Variation of Information: ', measuresAll['normalizedVariationInfo']
    # end if metric
    print '\n'
    
    if bShowEntropies:
        print 'entropy of first label set H(U)  = ', allEntropies['entropy1']
        print 'entropy of second label set H(V) = ', allEntropies['entropy2']
        print 'joint entropy H(U,V)             = ', allEntropies['joint']
        print 'conditional entropy H(U|V)       = ', allEntropies['conditional1']
        print 'conditional entropy H(V|U)       = ', allEntropies['conditional2']
        
# end main()

if __name__ == "__main__":
    # set the metrics to run 
    useMain = True # default is True to use command line interface; False to run standalone
    if useMain:
        main(sys.argv[1:])
         ######## example command lines: ########
        # 1. Normalized mutual information (NMI)
        # python quality_measures.py -d data_test/metrics/news20 -f news20_assignments_32_run1.csv -g news20_assignments_32_run1.csv -m NMI
        # output:
            # metric:  NMI
            # fileName1:  news20_assignments_32_run1.csv
            # fileName2:  news20_assignments_32_run1.csv
            # measuresAll:  {'normalizedMutualInfo': 1.0}
    
        # 2. Run all metrics (runall)
        # python quality_measures.py -d data_test/metrics/news20 -f news20_assignments_32_run1.csv -g news20_assignments_32_run1.csv -m all    
        # output:
            # metric:  runall
            # fileName1:  news20_assignments_32_run1.csv
            # fileName2:  news20_assignments_32_run4.csv
            # measuresAll:  {'variationInformation': 0.93453998956386597, 'informationDistance': 0.51534469787292103, 'normalizedVariationInfo': 0.25578066836602509, 'normalizedInfoDistance': 0.15932836127134131, 'normalizedMutualInfo': 0.85345235676699571, 'mutualInformation': 2.7191371844535559}
        
    elif not useMain: # run standalone
        # 3. test_measures unbound function can be run standalone without the command
        #       line interface; set useMain = False above
        clusterLabels1 = [0,1,3,0,6,5,4,14,10,23,22,11,17,20]
        clusterLabels2 = [1,0,5,4,3,0,9,19,23,26,21,19,20,28]
        #dAllMeas= test_measures(False, None, None, None, clusterLabels1, clusterLabels2)

        cMeas = measures(clusterLabels1,clusterLabels1)
        dAllMeas = cMeas.getAllMeasures()
        entPtr = cMeas.getEntropyPtr()
        allEntropies = entPtr.get_entropies()
        
        print 'entropy of first label set H(U)  = ', allEntropies['entropy1']
        print 'entropy of second label set H(V) = ', allEntropies['entropy2']
        print 'joint entropy H(U,V)             = ', allEntropies['joint']
        print 'conditional entropy H(U|V)       = ', allEntropies['conditional1']
        print 'conditional entropy H(V|U)       = ', allEntropies['conditional2']
        
        print 'mutualInformation        = ', dAllMeas['MI']
        print 'normalizedMutualInfo     = ', dAllMeas['NMI']
        print 'informationDistance      = ', dAllMeas['ID']
        print 'normalizedInfoDistance   = ', dAllMeas['NID']
        print 'variationInformation     = ', dAllMeas['VI']
        print 'normalizedVariationInfo1 = ', dAllMeas['NVI']
    else:
        print "check useMain setting\n"
        
            
# end if __name__
   
