# -*- coding: UTF-8 -*-
import os
def loadDataSet(filename):
    '''整理原始数据
       将原始数据整理成一个包含字典的列表
       格式是[ {'序列1'：'序列1类型'},
               {'序列2'：'序列2类型'},
               ...
               {'序列N'：'序列N类型'}]       
    '''    
    BASE_DIR = os.path.dirname(__file__) #获取当前文件夹的绝对路径
    file_path = os.path.join(BASE_DIR, filename)    
    f = open(file_path,"r")  
    data = []
    print "Load data..."    
    for line in f:
        lineSplit = line.split(",")
        gene = lineSplit[2].strip()
        dic = {gene:lineSplit[0]}
        data.append(dic)
    print "Load finished"    
    return data
    
def data2frozSet(data):
    '''整理第一步产生的字典列表，将字典的键值转变为frozenset
       将数据集转换为下面的格式
       例如[{frozenset(['0C','1C',...,'59G']):'EI'},
            {frozenset(['0A','1G',...,'59C']):'IE'}]
    '''
    frozDataList = []
    for item in data:
        gene = item.keys()[0]
        geneType = item.values()[0]
        genefSet = frozenset([])
        for count in range(len(gene)):
            genePoisition = str(count) + gene[count]
            tempfSet = frozenset([genePoisition])
            genefSet = genefSet.union(tempfSet)
        dic = {genefSet: geneType}
        frozDataList.append(dic)
    return frozDataList

def creatC1(frozData):
    '''产生 1-频繁项候选集
       输入参数是上一步产生的frozData
	 结果是一个列表，包含所有不重复的1项集
       例如[frozenset(['0C']), frozenset(['1C']) ...frozenset(['59T'])]
    '''
    C1 = []
    for data in frozData:
        gene = data.keys()[0]
        for item in gene:
            if not [item] in C1:
                C1.append([item])             
    C1.sort()
    return map(frozenset, C1)
    
def Ck2Lk(froData, candidate, minSupport = 0.1):
    '''从Ck(k-频繁项候选集)生产Lk(k-频繁项集)
       输入参数是原始的frozData，候选集candidate以及minSupport，默认为0.1
       输出为Lk,同时输出一个包含支持度的字典Cksupp
    '''
    ssCnt = {}  
    for data in froData:
        frozData = data.keys()[0]
        for can in candidate:
            if can.issubset(frozData):
                if not ssCnt.has_key(can):
                    ssCnt[can] = 1
                else:
                    ssCnt[can] += 1
    dataSize = float(len(frozData))   
    Lk = []
    CkSupp = {}
    print "Computing support..."
    for key in ssCnt:
        support = ssCnt[key]/dataSize
        if support >= minSupport:
            Lk.insert(0, key)
            #Lk.append(key)
        CkSupp[key] = support
    return Lk, CkSupp
                        
def TBARgenCk(Lk, L1):
    '''从Lk(k-频繁项集)产生C(k+1)(k+1频繁项候选集)
       输入为Lk和L1，输出为CkplusList
    '''
    LkSize = len(Lk)
    L1Size = len(L1)
    CkplusList = []
    for count1 in range(LkSize):
        FroGenePos1 = Lk[count1]
        GenePosList = list(FroGenePos1)
        maxPos = -1
        for GenePos in GenePosList:
            Pos = int(GenePos[:-1]) #找出字符串最后一位之前的所有位
            if Pos > maxPos:
                maxPos = Pos
        for count2 in range(L1Size):
            FroGenePos2 = L1[count2]
            GenePos2 = list(FroGenePos2)[0] #L1为1项集，只取第一个元素，也是最后一个元素
            Pos2 = int(GenePos2[:-1])
            if maxPos < Pos2:
                Ckplus = FroGenePos1.union(FroGenePos2)
                CkplusList.append(Ckplus)
    return CkplusList
    
def genRule(Lk,froData):
    '''计算置信度产生规则，每次取计算出来的最大的置信度减支持度作为置信度
       输入为Lk和最初的froData
       输出为rule和confDic,rule为一个包含支持度的字典,confDic为一个字典，key为frozenset的rule
    '''
    confDic = {}
    NCount = 0.0
    EICount = 0.0
    IECount = 0.0
    print "Match Lk with froData"
    for data in froData:
        geneType = data.values()[0]
        if geneType == 'N':
            NCount += 1
        elif geneType == 'EI':
            EICount += 1
        elif geneType == 'IE':
            IECount += 1
            
    maxConf = 0   
    for item in Lk:               
        Allsupport = 0.0
        Nsupport = 0.0
        EIsupport = 0.0
        IEsupport = 0.0 
        for data in froData:
            gene = data.keys()[0]
            geneType = data.values()[0]
            if item.issubset(gene):
                Allsupport += 1
                if geneType == 'N':
                    Nsupport += 1
                elif geneType == 'EI':
                    EIsupport += 1
                elif geneType == 'IE':
                    IEsupport += 1
        
        confDic[item] = {'N'     : Nsupport/(NCount+1),
                         'EI'    : EIsupport/(EICount+1),
                         'IE'    : IEsupport/(IECount+1),
                         'total' : Allsupport}       
        maxConf = max(maxConf, confDic[item]['N'], confDic[item]['EI'], confDic[item]['IE'])
    minConf = maxConf - 0.05
    rule = []
    #print 'Selecting rules...'
    for key in confDic:
        maxConfselect = max(confDic[key]['N'], confDic[key]['EI'], confDic[key]['IE'])
        if maxConfselect >= minConf:
            if confDic[key]['N'] == maxConfselect:
                rule.append({key: {'N': maxConfselect}})
            elif confDic[key]['EI'] == maxConfselect:
                rule.append({key: {'EI': maxConfselect}})
            elif confDic[key]['IE'] == maxConfselect:
                rule.append({key: {'IE': maxConfselect}})
    return rule, confDic