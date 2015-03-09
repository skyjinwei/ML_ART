# -*- coding: utf-8 -*-
import ARTNewJ
data = ARTNewJ.loadDataSet("spliceTrain1.txt")
frozData = ARTNewJ.data2frozSet(data)
C = []  #候选频繁项集全局变量
L = []  #频繁项集列表全部变量
CSupp = []  #带支持度的候选频繁项集全部变量
rule = []   #规则列表全部变量
confDic = []  #confDic为一个字典，key为rule
calculate = [] #计算第i次迭代所覆盖的数据

'''生成规则，论文中最多只有2项集的规则，在这里最多生成3项集
   输入参数为未被覆盖的frozData，改变全局变量C,L,rule等
'''
def ruleGen(frozData):
    C1 = ARTNewJ.creatC1(frozData)
    L1, C1Supp = ARTNewJ.Ck2Lk(frozData, C1)
    rule1, confDic1 = ARTNewJ.genRule(L1, frozData)
    C.append(C1); L.append(L1); CSupp.append(C1Supp)
    rule.append(rule1); confDic.append(confDic1)
    count = 1
    while(count < 2 and len(C[count - 1]) > 0):
        print "Selecting rules with "+str(count+1)+" sets"    
        Ctemp = ARTNewJ.TBARgenCk(L[count - 1],L1)
        C.append(Ctemp)
        Ltemp, Csupptemp = ARTNewJ.Ck2Lk(frozData, C[count])
        L.append(Ltemp)
        CSupp.append(Csupptemp)
        ruletemp, confDictemp = ARTNewJ.genRule(L[count], frozData)
        rule.append(ruletemp)
        confDic.append(confDictemp)
        count += 1
    
'''对输入的数据进行规则的遍历，并返回未被当前规则覆盖的数据
   输入是frozendata,输出也是frozendata
'''
def ruleTraversal(frozData, number):
    uncoverData = []
    count = 1;
    for data in frozData:
        if(judge(data)):
            count += 1
        else:
            uncoverData.append(data)
    calculate.append({number:count})
    return uncoverData
    
'''判断某条数据是否在当前规则中
   在则返回true，不在则返回false
'''  
def judge(data):
    for rulek in rule:
        for ruleitem in rulek:
            if ruleitem.keys()[0].issubset(data.keys()[0]):
                return True
    return False
                    
'''输入测试数据集，把计算出来的误分率加到原规则中去
'''
def test():
    dataTest = ARTNewJ.loadDataSet("spliceTest.txt")
    frozDataTest = ARTNewJ.data2frozSet(dataTest)
    count = 0.0
    for dataT in frozDataTest:
        count += 1
    print "count:", count
    for rulek in rule:
        for ruleitem in rulek:
            err = 0.0
            for dataT in frozDataTest:
                if ruleitem.keys()[0].issubset(dataT.keys()[0]):
                    if ruleitem.values()[0].keys()[0] != dataT.values()[0]: #如果规则的type和数据的type不相同，则预测错误
                        err += 1
            ruleitem.values()[0]['err'] = err/count  #在原来规则的基础上加上误分率

'''main函数初始化将原始frozData赋值给局部变量uncoverfrozData
   如果uncoverfrozData非空，则继续生成规则，并移除数据
'''        
def main():
    uncoverfrozData = frozData
    number = 0
    while(len(uncoverfrozData)):
        number += 1
        print "The " + str(number) + "th traversal"
        ruleGen(uncoverfrozData)
        uncoverfrozData = ruleTraversal(uncoverfrozData, number)
    test()
           
if __name__ == '__main__':
    main()
    
