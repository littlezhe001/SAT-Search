# -*- coding: utf-8 -*-


import os
import time
import random

FullRound = 32

SearchRoundStart = 1
SearchRoundEnd = 32
InitialLowerBound = 0

GroupConstraintChoice = 1

# Parameters for choice 1
GroupNumForChoice1 = 1

DiffActiveSbox = list([])
for i in range(FullRound):
    DiffActiveSbox += [0]


def CountClausesInRoundFunction(Round, ActiveSbox, clause_num):  # 计算在轮函数中的子句数量
    count = clause_num
    # Nonzero input
    count += 1
    # Cluases for Sbox
    for r in range(Round):
        for i in range(16):
            for j in range(43):
                count += 1
    return count


def CountClausesInSequentialEncoding(main_var_num, cardinalitycons, clause_num):  # 计算子句在序列编码中的数量   对应着11页的公式
    count = clause_num
    n = main_var_num  # S盒的总数
    k = cardinalitycons  # 活跃S盒的上界
    if (k > 0):
        count += 1
        for j in range(1, k):
            count += 1
        for i in range(1, n - 1):
            count += 3
        for j in range(1, k):
            for i in range(1, n - 1):
                count += 2
        count += 1
    if (k == 0):
        for i in range(n):
            count += 1
    return count


# 产生序列编码
def GenSequentialEncoding(x, u, main_var_num, cardinalitycons,
                          fout):  # x代表11页公式（1）中的x即S盒的数量,u代表s即引入的辅助变量，cardinalitycons表示所设置的活跃S盒最小数量
    n = main_var_num  # 为总S盒数量
    k = cardinalitycons
    # print("x:  "+str(x))
    if (k > 0):
        clauseseq = "-" + str(x[0] + 1) + " " + str(u[0][0] + 1) + " 0" + "\n"  # 第一行
        fout.write(clauseseq)
        for j in range(1, k):  # 第二行
            clauseseq = "-" + str(u[0][j] + 1) + " 0" + "\n"
            fout.write(clauseseq)
        for i in range(1, n - 1):
            clauseseq = "-" + str(x[i] + 1) + " " + str(u[i][0] + 1) + " 0" + "\n"  # 第三行
            fout.write(clauseseq)
            clauseseq = "-" + str(u[i - 1][0] + 1) + " " + str(u[i][0] + 1) + " 0" + "\n"  # 第四行
            fout.write(clauseseq)
            clauseseq = "-" + str(x[i] + 1) + " " + "-" + str(u[i - 1][k - 1] + 1) + " 0" + "\n"  # 第七行
            fout.write(clauseseq)
        for j in range(1, k):  # 第五行和第六行
            for i in range(1, n - 1):
                clauseseq = "-" + str(x[i] + 1) + " " + "-" + str(u[i - 1][j - 1] + 1) + " " + str(
                    u[i][j] + 1) + " 0" + "\n"
                fout.write(clauseseq)
                clauseseq = "-" + str(u[i - 1][j] + 1) + " " + str(u[i][j] + 1) + " 0" + "\n"
                fout.write(clauseseq)

        clauseseq = "-" + str(x[n - 1] + 1) + " " + "-" + str(u[n - 2][k - 1] + 1) + " 0" + "\n"  # 最后一行
        fout.write(clauseseq)
    if (k == 0):
        for i in range(n):
            clauseseq = "-" + str(x[i] + 1) + " 0" + "\n"
            fout.write(clauseseq)


def CountClausesForMatsuiStrategy(n, k, left, right, m, clausenum):  # left相当与14页的e1,right相当于e2
    count = clausenum
    if (m > 0):
        if ((left == 0) and (right < n - 1)):
            for i in range(1, right + 1):  # 公式（7）
                count += 1
        if ((left > 0) and (right == n - 1)):  # 公式（9）
            for i in range(0, k - m):
                count += 1
            for i in range(0, k - m + 1):
                count += 1
        if ((left > 0) and (right < n - 1)):  # 公式（8）
            for i in range(0, k - m):
                count += 1
    if (m == 0):
        for i in range(left, right + 1):
            count += 1
    return count


def GenMatsuiConstraint(x, u, n, k, left, right, m, fout):  # 产生约束句子
    # print("m is :"+str(m)+"  left :" + str(left) + "  right:" + str(right) + " x length:"+ str(len(x)))
    if (m > 0):
        if ((left == 0) and (right < n - 1)):
            for i in range(1, right + 1):
                clauseseq = "-" + str(x[i] + 1) + " " + "-" + str(u[i - 1][m - 1] + 1) + " 0" + "\n"
                fout.write(clauseseq)
        if ((left > 0) and (right == n - 1)):
            for i in range(0, k - m):
                clauseseq = str(u[left - 1][i] + 1) + " " + "-" + str(u[right - 1][i + m] + 1) + " 0" + "\n"
                fout.write(clauseseq)
            for i in range(0, k - m + 1):
                clauseseq = str(u[left - 1][i] + 1) + " " + "-" + str(x[right] + 1) + " " + "-" + str(
                    u[right - 1][i + m - 1] + 1) + " 0" + "\n"
                fout.write(clauseseq)
        if ((left > 0) and (right < n - 1)):
            for i in range(0, k - m):
                clauseseq = str(u[left - 1][i] + 1) + " " + "-" + str(u[right][i + m] + 1) + " 0" + "\n"
                fout.write(clauseseq)
    if (m == 0):
        for i in range(left, right + 1):
            clauseseq = "-" + str(x[i] + 1) + " 0" + "\n"
            fout.write(clauseseq)
def CountGenSboxConstrain(w,Round,clausenum):
    count = clausenum
    for r in range(1,Round):
        for i in range(16):
            count+=1
    return count

def GenSboxConstrain(w,Round,fout):
    p = [0, 4, 8, 12, 0, 4, 8, 12, 0, 4, 8, 12, 0, 4, 8, 12, 1, 5, 9, 13, 1, 5, 9, 13, 1, 5, 9, 13, 1, 5, 9, 13, 2, 6, 10,
     14, 2, 6, 10, 14, 2, 6, 10, 14, 2, 6, 10, 14, 3, 7, 11, 15, 3, 7, 11, 15, 3, 7, 11, 15, 3, 7, 11, 15]
    for r in range(1,Round):
        for i in range(16):
            X = []
            for j in range(4):
                X.append(p[i*4+j])
            X.append(w[r][i])
            clauseseq = str(X[0] + 1) + " " + str(X[1] + 1) + " " + str(X[2] + 1) + " " + str(X[3] + 1) + " "+ "-" + str(X[4] + 1) + " 0" + "\n"
            fout.write(clauseseq)

def Decision(Round, ActiveSbox, MatsuiRoundIndex, MatsuiCount, flag):
    TotalSbox = 16 * Round
    count_var_num = 0  # 变量的数量
    # Declare variables
    xin = []  # 进入S盒的比特  前面有-，代表该比特为0，
    w = []  # 存储每一轮的活跃S盒数量
    xout = []  # 输出的x
    for i in range(Round):
        xin.append([])
        w.append([])
        xout.append([])
        for j in range(64):
            xin[i].append(0)
        for j in range(16):
            w[i].append(0)
        for j in range(64):
            xout[i].append(0)
    # Allocate variables
    for i in range(Round):
        for j in range(64):
            xin[i][j] = count_var_num
            count_var_num += 1
    for i in range(Round - 1):
        for j in range(64):
            xout[i][j] = xin[i + 1][j]
    for i in range(64):
        xout[Round - 1][i] = count_var_num
        count_var_num += 1
    for i in range(Round):
        for j in range(16):
            w[i][j] = count_var_num
            count_var_num += 1

    auxiliary_var_u = []  # 辅助变量
    for i in range(TotalSbox - 1):
        auxiliary_var_u.append([])
        for j in range(ActiveSbox):
            auxiliary_var_u[i].append(count_var_num)  # 用于序列编码那里
            count_var_num += 1

    # Count the number of clauses in the round function
    count_clause_num = 0
    count_clause_num = CountClausesInRoundFunction(Round, ActiveSbox, count_clause_num)
    # Count the number of clauses in the original sequential encoding
    Main_Var_Num = 16 * Round
    CardinalityCons = ActiveSbox
    count_clause_num = CountClausesInSequentialEncoding(Main_Var_Num, CardinalityCons, count_clause_num)  # 计算序列编码的句子
    # Count the number of clauses for Matsui's strategy
    # for matsui_count in range(0, MatsuiCount):
    #     StartingRound = MatsuiRoundIndex[matsui_count][0]
    #     EndingRound = MatsuiRoundIndex[matsui_count][1]
    #     LeftNode = 16 * StartingRound
    #     RightNode = 16 * EndingRound - 1
    #     PartialCardinalityCons = ActiveSbox - DiffActiveSbox[StartingRound] - DiffActiveSbox[Round - EndingRound]
    #     count_clause_num = CountClausesForMatsuiStrategy(Main_Var_Num, CardinalityCons, LeftNode, RightNode,
    #                                                      PartialCardinalityCons, count_clause_num)
    # count_clause_num = CountGenSboxConstrain(w,Round,count_clause_num)
    # Open file
    file = open("Problem-Round" + str(Round) + "-Active" + str(ActiveSbox) + ".cnf", "w")
    print("所需变量个数为：" + str(count_var_num) + " " + "所需要的句子数量为：" + str(count_clause_num))
    file.write("p cnf " + str(count_var_num) + " " + str(count_clause_num) + "\n")  # 第一行写变量的总数   句子的数量
    # Add constraints to claim nonzero input difference
    clauseseq = ""
    for i in range(64):
        clauseseq += str(xin[0][i] + 1) + " "
    clauseseq += "0" + "\n"
    file.write(clauseseq)

    # Add constraints for the round function   对轮函数进行约束
    for r in range(Round):
        y = list([])
        P = [0, 16, 32, 48, 1, 17, 33, 49, 2, 18, 34, 50, 3, 19, 35, 51, 4, 20, 36, 52, 5, 21, 37, 53, 6, 22, 38, 54, 7,
             23, 39, 55, 8, 24, 40, 56, 9, 25, 41, 57, 10, 26, 42, 58, 11, 27, 43, 59, 12, 28, 44, 60, 13, 29, 45, 61,
             14, 30, 46, 62, 15, 31, 47, 63]
        SymbolicCNFConstraintForSbox = [  # Differential PRESENT (43)
            [1, 0, 1, 9, 0, 1, 1, 9, 9], [9, 0, 0, 9, 1, 0, 1, 0, 9], [0, 9, 9, 1, 1, 1, 1, 0, 9],
            [9, 1, 0, 0, 9, 0, 1, 1, 9], [1, 0, 1, 9, 1, 1, 0, 0, 9], [1, 1, 0, 9, 0, 1, 1, 0, 9],
            [1, 1, 0, 9, 1, 1, 0, 9, 9], [0, 0, 1, 1, 1, 9, 0, 9, 9], [9, 1, 0, 1, 1, 1, 1, 9, 9],
            [9, 0, 1, 0, 1, 9, 1, 1, 9], [9, 1, 0, 0, 1, 1, 9, 1, 9], [0, 0, 0, 9, 1, 9, 1, 9, 9],
            [9, 0, 1, 0, 9, 0, 0, 1, 9], [9, 1, 0, 0, 0, 0, 9, 1, 9], [9, 1, 1, 0, 0, 9, 0, 9, 9],
            [9, 0, 1, 0, 9, 1, 1, 1, 9], [0, 0, 1, 9, 1, 0, 0, 9, 9], [9, 1, 1, 1, 1, 0, 1, 9, 9],
            [9, 0, 1, 1, 1, 1, 1, 9, 9], [0, 9, 0, 9, 9, 0, 0, 0, 1], [0, 1, 0, 9, 0, 0, 1, 9, 9],
            [0, 1, 1, 9, 0, 9, 0, 0, 9], [1, 1, 9, 0, 1, 9, 1, 1, 9], [0, 1, 0, 1, 9, 1, 1, 9, 9],
            [0, 1, 1, 9, 0, 9, 1, 1, 9], [9, 1, 1, 0, 1, 9, 1, 0, 9], [0, 1, 1, 9, 1, 9, 0, 1, 9],
            [9, 0, 0, 0, 1, 9, 9, 0, 9], [9, 9, 9, 9, 0, 0, 0, 0, 1], [9, 9, 9, 1, 0, 1, 0, 1, 9],
            [1, 9, 1, 1, 0, 9, 1, 9, 9], [0, 0, 9, 9, 0, 0, 9, 0, 1], [1, 0, 0, 1, 9, 9, 9, 1, 9],
            [1, 1, 9, 1, 1, 9, 0, 9, 9], [9, 0, 0, 9, 0, 9, 0, 1, 9], [9, 9, 9, 0, 0, 1, 0, 0, 9],
            [9, 1, 9, 9, 9, 9, 9, 9, 0], [0, 0, 0, 0, 9, 9, 9, 9, 1], [0, 0, 0, 1, 9, 9, 9, 0, 9],
            [9, 0, 0, 0, 9, 9, 1, 0, 9], [9, 9, 1, 9, 9, 9, 9, 9, 0], [9, 9, 9, 9, 9, 9, 9, 1, 0],
            [1, 9, 9, 9, 9, 9, 9, 9, 0]]
        for i in range(64):
            y += [xout[r][P[i]]]  # 相当于把下一轮的x进行p置换   感觉应该进行P逆
        for i in range(16):
            for j in range(43):
                X = list([])
                for k in range(4):
                    X += [xin[r][4 * i + k]]
                for k in range(4):
                    X += [y[4 * i + k]]
                X += [w[r][i]]
                clauseseq = ""
                for k in range(9):
                    if (SymbolicCNFConstraintForSbox[j][k] == 1):
                        clauseseq += "-" + str(X[k] + 1) + " "
                    if (SymbolicCNFConstraintForSbox[j][k] == 0):
                        clauseseq += str(X[k] + 1) + " "
                clauseseq += "0" + "\n"
                file.write(clauseseq)
    # Add constraints for the original sequential encoding
    Main_Vars = list([])
    for r in range(Round):
        for i in range(16):
            Main_Vars += [w[r][i]]
    # print(w)
    GenSequentialEncoding(Main_Vars, auxiliary_var_u, Main_Var_Num, CardinalityCons,
                          file)  # Main_Vars相当于x  auxiliary_var_u相当于s

    # GenSboxConstrain(w,Round,file)
    # Add constraints for Matsui's strategy
    # for matsui_count in range(0, MatsuiCount):
    #     StartingRound = MatsuiRoundIndex[matsui_count][0]
    #     EndingRound = MatsuiRoundIndex[matsui_count][1]
    #     LeftNode = 16 * StartingRound
    #     RightNode = 16 * EndingRound - 1
    #     PartialCardinalityCons = ActiveSbox - DiffActiveSbox[StartingRound] - DiffActiveSbox[Round - EndingRound]
    #     GenMatsuiConstraint(Main_Vars, auxiliary_var_u, Main_Var_Num, CardinalityCons, LeftNode, RightNode,
    #                         PartialCardinalityCons, file)

    file.close()
    time_start = time.time()
    # Call solver cadical
    order = "cadical " + "Problem-Round" + str(Round) + "-Active" + str(ActiveSbox) + ".cnf > Round" + str(
        Round) + "-Active" + str(ActiveSbox) + "-solution.out"
    os.system(order)
    # Extracting results
    order = "sed -n '/s SATISFIABLE/p' Round" + str(Round) + "-Active" + str(
        ActiveSbox) + "-solution.out > SatSolution.out"
    os.system(order)
    order = "sed -n '/s UNSATISFIABLE/p' Round" + str(Round) + "-Active" + str(
        ActiveSbox) + "-solution.out > UnsatSolution.out"
    os.system(order)
    satsol = open("SatSolution.out")
    unsatsol = open("UnsatSolution.out")
    satresult = satsol.readlines()
    unsatresult = unsatsol.readlines()
    satsol.close()
    unsatsol.close()
    if ((len(satresult) == 0) and (len(unsatresult) > 0)):
        flag = False
    if ((len(satresult) > 0) and (len(unsatresult) == 0)):
        flag = True
    order = "rm SatSolution.out"
    os.system(order)
    order = "rm UnsatSolution.out"
    os.system(order)
    # Removing cnf file
    order = "rm Problem-Round" + str(Round) + "-Active" + str(ActiveSbox) + ".cnf"
    os.system(order)
    time_end = time.time()
    # Printing solutions
    if (flag == True):
        print(
            "Round:" + str(Round) + "; Active: " + str(ActiveSbox) + "; Sat; TotalCost: " + str(time_end - time_start))
    else:
        print("Round:" + str(Round) + "; Active: " + str(ActiveSbox) + "; Unsat; TotalCost: " + str(
            time_end - time_start))
    return flag


# main function
# 为了自动计算最优开始约束轮数
maxd = 0
CountSbox = InitialLowerBound  # 初始预估活跃S盒上界
TotalTimeStart = time.time()
for totalround in range(SearchRoundStart, SearchRoundEnd):
    flag = False
    time_start = time.time()
    MatsuiRoundIndex = []
    MatsuiCount = 0
    # Generate Matsui condition under choice 1
    if (GroupConstraintChoice == 1):
        for group in range(0, GroupNumForChoice1):
            for round in range(1, totalround - group - 1):
                MatsuiRoundIndex.append([])
                MatsuiRoundIndex[MatsuiCount].append(group)
                MatsuiRoundIndex[MatsuiCount].append(group + round)
                MatsuiCount += 1
    # print(MatsuiRoundIndex)
    # Printing Matsui conditions
    file = open("MatsuiCondition.out", "a")
    resultseq = "Round: " + str(totalround) + "; Partial Constraint Num: " + str(MatsuiCount) + "\n"
    file.write(resultseq)
    file.write(str(MatsuiRoundIndex) + "\n")
    file.close()


    while (flag == False):
        flag = Decision(totalround, CountSbox, MatsuiRoundIndex, MatsuiCount, flag)
        CountSbox += 1
    DiffActiveSbox[totalround] = CountSbox - 1

    time_end = time.time()
    file = open("RunTimeSummarise.out", "a")
    resultseq = "Round: " + str(totalround) + "; Active S-box: " + str(
        DiffActiveSbox[totalround]) + "; Runtime: " + str(time_end - time_start) + "\n"
    file.write(resultseq)
    file.close()
    print(time_end - TotalTimeStart)
print(str(DiffActiveSbox))
TotalTimeEnd = time.time()
print("Total Runtime: " + str(TotalTimeEnd - TotalTimeStart))
file = open("RunTimeSummarise.out", "a")
resultseq = "Total Runtime: " + str(TotalTimeEnd - TotalTimeStart)
file.write(resultseq)
