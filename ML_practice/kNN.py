from numpy import *
import operator

def creatDataSet():
    group = array([[1.0, 1.1], [1.0, 1.0], [0,0],[0, 0.1]])
    labels = ['A', 'A', 'B', 'B']
    return group, labels

group, labels = creatDataSet()
print group

def classify0(inX, dataSet, labels, k):

    dataSetSize = dataSet.shape[0]
    print dataSetSize
    diffMat = tile(inX, (dataSetSize, 1)) - dataSet
    sqDiffMat = diffMat**2
    sqDistances = sqDiffMat.sum(axis=1)
    distances = sqDistances**0.5
    print dataSetSize, distances
    sortedDistIndicies = distances.argsort()
    print sortedDistIndicies, labels
    classCount = {}
    for i in range(k):
        print i
        voteIlabel = labels[sortedDistIndicies[i]]
        print voteIlabel, classCount.get(voteIlabel, 0)
        classCount[voteIlabel] = classCount.get(voteIlabel, 0) + 1
        print classCount
    sortedClassCount = sorted(classCount.iteritems(), key=operator.itemgetter(0), reverse=False)

    print sortedClassCount
    return sortedClassCount[0][0]
inX, k = [0, 0], 3
classify0(inX, group, labels, k)
