import os
from list import *
import time

FullRound = 43

SearchRoundStart = 1
SearchRoundEnd = 14
InitialLowerBound = 1
GroupConstraintChoice = 1
# Parameters for choice 10
GroupNumForChoice1 = 1

DiffActiveSbox = list([])
for i in range(FullRound):
    DiffActiveSbox += [0]


def CountClausesInRoundFunction(Round, ActiveSbox, clause_num):
    count = clause_num
    # Nonzero input
    count += 1
    # Cluases for Round function
    for r in range(Round):
        for block in range(16):
            count += 36
        for i in range(32):
            count += 40
    return count


# 新的序列编码句子计数
def CountClausesInSequentialEncoding(main_var_num, cardinalitycons, clause_num):  # 计算子句在序列编码中的数量   对应着11页的公式
    count = clause_num
    n = main_var_num  # S盒的总数
    k = cardinalitycons  # 活跃S盒的上界
    m = n - k
    if (k > 0):
        count += 1
        for i in range(1, n - 1):
            count += 2

        for i in range(k, n - 1):
            count += 1
        for i in range(1, n - 1):
            if i < k:
                for j in range(1, i):
                    count += 2
                count += 1
            else:
                for j in range(1, k):
                    count += 2
        count += 1
    if (k == 0):
        for i in range(n):
            count += 1
    return count


# 新的序列编码
def GenSequentialEncoding(x, u, main_var_num, cardinalitycons,
                          fout):  # x代表11页公式（1）中的x即S盒的数量,u代表s即引入的辅助变量，cardinalitycons表示所设置的活跃S盒最小数量
    n = main_var_num  # 为总S盒数量
    k = cardinalitycons
    if (k > 0):
        clauseseq = "-" + str(x[0] + 1) + " " + str(u[0][0] + 1) + " 0" + "\n"  # 第一行
        fout.write(clauseseq)

        for i in range(1, n - 1):
            clauseseq = "-" + str(x[i] + 1) + " " + str(u[i][0] + 1) + " 0" + "\n"  # 第三行
            fout.write(clauseseq)
            clauseseq = "-" + str(u[i - 1][0] + 1) + " " + str(u[i][0] + 1) + " 0" + "\n"  # 第四行
            fout.write(clauseseq)
        for i in range(k, n - 1):
            clauseseq = "-" + str(x[i] + 1) + " " + "-" + str(u[i - 1][k - 1] + 1) + " 0" + "\n"  # 第七行
            fout.write(clauseseq)
        for i in range(1, n - 1):
            if i < k:
                for j in range(1, i):
                    clauseseq = "-" + str(x[i] + 1) + " " + "-" + str(u[i - 1][j - 1] + 1) + " " + str(
                        u[i][j] + 1) + " 0" + "\n"
                    fout.write(clauseseq)
                    clauseseq = "-" + str(u[i - 1][j] + 1) + " " + str(u[i][j] + 1) + " 0" + "\n"
                    fout.write(clauseseq)
                clauseseq = "-" + str(x[i] + 1) + " " + "-" + str(u[i - 1][i - 1] + 1) + " " + str(
                    u[i][i] + 1) + " 0" + "\n"
                fout.write(clauseseq)
            else:
                for j in range(1, k):
                    clauseseq = "-" + str(x[i] + 1) + " " + "-" + str(u[i - 1][j - 1] + 1) + " " + str(
                        u[i][j] + 1) + " 0" + "\n"
                    fout.write(clauseseq)
                    clauseseq = "-" + str(u[i - 1][j] + 1) + " " + str(u[i][j] + 1) + " 0" + "\n"
                    fout.write(clauseseq)

        clauseseq = "-" + str(x[n - 1] + 1) + " " + "-" + str(u[n - 2][k - 1] + 1) + " 0" + "\n"  # 最后一行
        fout.write(clauseseq)


def CountClausesForMatsuiStrategy(n, k, left, right, m, clausenum):
    count = clausenum
    if (m > 0):
        if ((left == 0) and (right < n - 1)):
            for i in range(m, right + 1):
                count += 1
        if ((left > 0) and (right == n - 1)):
            for i in range(0, k - m):
                count += 1
            for i in range(0, k - m + 1):
                count += 1
        if ((left > 0) and (right < n - 1)):
            for i in range(0, k - m):
                count += 1
    if (m == 0):
        for i in range(left, right + 1):
            count += 1
    return count


def GenMatsuiConstraint(x, u, n, k, left, right, m, fout):
    if (m > 0):
        if ((left == 0) and (right < n - 1)):
            for i in range(m, right + 1):
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


def genclause(X, i, l1,flag=False):
    ibin = bin(i)[2:]
    while len(ibin) != l1:
        ibin = '0' + ibin
    clause = ''
    for i in range(l1):
        if ibin[i] == '0':
            clause = clause + str(X[i] + 1) + ' '
        else:
            clause = clause + '-' + str(X[i] + 1) + ' '
    if flag:
        w = 1
        for i in ibin:
            w = w ^ int(i, 2)
        if w:
            clause = clause + '-' + str(X[-1] + 1) + ' '
        else:
            clause = clause + str(X[-1] + 1) + ' '
    clause = clause + '0' + '\n'
    return clause


def Decision(Round, ActiveSbox, MatsuiRoundIndex, MatsuiCount, flag):
    TotalSbox = 16 * Round
    count_var_num = 0
    time_start = time.time()
    xin = []
    xout = []
    asout = []
    dsout = []
    w = []
    for i in range(Round):
        xin.append([])
        xout.append([])
        asout.append([])
        dsout.append([])
        w.append([])
        for j in range(128):
            xin[i].append(0)
            xout[i].append(0)
        for j in range(32):
            asout[i].append(0)
            dsout[i].append(0)
        for j in range(16):
            w[i].append(0)

    # Allocate variables
    for i in range(Round):
        for j in range(128):
            xin[i][j] = count_var_num
            count_var_num += 1
        for j in range(32):
            asout[i][j] = count_var_num
            count_var_num += 1
        for j in range(32):
            dsout[i][j] = count_var_num
            count_var_num += 1
        for j in range(16):
            w[i][j] = count_var_num
            count_var_num += 1
    for i in range(Round - 1):
        for j in range(128):
            xout[i][j] = xin[i + 1][j]
    for i in range(128):
        xout[Round - 1][i] = count_var_num
        count_var_num += 1
    auxiliary_var_u = []
    for i in range(TotalSbox - 1):
        auxiliary_var_u.append([])
        for j in range(ActiveSbox):
            auxiliary_var_u[i].append(count_var_num)
            count_var_num += 1

    # Count the number of clauses in the round function
    count_clause_num = 0
    count_clause_num = CountClausesInRoundFunction(Round, ActiveSbox, count_clause_num)
    # Count the number of clauses in the original sequential encoding
    Main_Var_Num = 16 * Round
    CardinalityCons = ActiveSbox
    count_clause_num = CountClausesInSequentialEncoding(Main_Var_Num, CardinalityCons, count_clause_num)
    # Count the number of clauses for Matsui's strategy
    for matsui_count in range(0, MatsuiCount):
        StartingRound = MatsuiRoundIndex[matsui_count][0]
        EndingRound = MatsuiRoundIndex[matsui_count][1]
        LeftNode = 16 * StartingRound
        RightNode = 16 * EndingRound - 1
        PartialCardinalityCons = ActiveSbox - DiffActiveSbox[StartingRound] - DiffActiveSbox[Round - EndingRound]
        count_clause_num = CountClausesForMatsuiStrategy(Main_Var_Num, CardinalityCons, LeftNode, RightNode,
                                                         PartialCardinalityCons, count_clause_num)

    # Open file
    file = open("Problem-Round" + str(Round) + "-Active" + str(ActiveSbox) + ".cnf", "w")
    file.write("p cnf " + str(count_var_num) + " " + str(count_clause_num) + "\n")
    # Add constraints to claim nonzero input difference
    clauseseq = ""
    for i in range(128):
        clauseseq += str(xin[0][i] + 1) + " "
    clauseseq += "0" + "\n"
    file.write(clauseseq)

    # Add constraints for the round function

    for r in range(Round):
        # 限制S盒
        for block in range(8):
            X = []
            for i in range(4):
                X += [xin[r][i * 8 + block]]
            for i in range(4):
                X += [asout[r][i * 8 + block]]
            X += [w[r][block]]
            for i in range(36):
                clauseseq = ""
                for k in range(9):
                    if SymbolicCNFConstraintForSbox[i][k] == 1:
                        clauseseq += "-" + str(X[k] + 1) + " "
                    if SymbolicCNFConstraintForSbox[i][k] == 0:
                        clauseseq += str(X[k] + 1) + " "
                clauseseq += "0" + "\n"
                file.write(clauseseq)
        for block in range(8):
            X = []
            for i in range(4):
                X += [xin[r][i * 8 + block + 96]]
            for i in range(4):
                X += [dsout[r][i * 8 + block]]
            X += [w[r][block + 8]]
            for i in range(36):
                clauseseq = ""
                for k in range(9):
                    if SymbolicCNFConstraintForSbox[i][k] == 1:
                        clauseseq += "-" + str(X[k] + 1) + " "
                    if SymbolicCNFConstraintForSbox[i][k] == 0:
                        clauseseq += str(X[k] + 1) + " "
                clauseseq += "0" + "\n"
                file.write(clauseseq)
        # 轮函数线性操作  3  10
        xouta = xout[r][32:64]
        xoutb = xout[r][0:32]
        xoutc = xout[r][96:128]
        xoutd = xout[r][64:96]
        for k in range(2):
            if k == 0:
                xorinput1 = asout[r][3:] + asout[r][:3]
                xorinput2 = asout[r][10:] + asout[r][:10]
                xorinput3 = asout[r]
                xorinput4 = xin[r][32:64]
                xorinput5 = xoutb

            else:
                xorinput1 = dsout[r][3:] + dsout[r][:3]
                xorinput2 = dsout[r][10:] + dsout[r][:10]
                xorinput3 = dsout[r]
                xorinput4 = xin[r][64:96]
                xorinput5 = xoutc
            for i in range(32):
                a = xorinput1[i]
                b = xorinput2[i]
                c = xorinput3[i]
                d = xorinput4[i]
                e = xorinput5[i]
                clauseseq = str(a + 1) + " " + str(b + 1) + " " + str(c + 1) + " " + str(d + 1) + " " + "-" + str(
                    e + 1) + " " + "0" + "\n"
                file.write(clauseseq)
                clauseseq = str(a + 1) + " " + str(b + 1) + " " + str(c + 1) + " " + "-" + str(d + 1) + " " + str(
                    e + 1) + " " + "0" + "\n"
                file.write(clauseseq)
                clauseseq = str(a + 1) + " " + str(b + 1) + " " + "-" + str(c + 1) + " " + str(d + 1) + " " + str(
                    e + 1) + " " + "0" + "\n"
                file.write(clauseseq)
                clauseseq = str(a + 1) + " " + "-" + str(b + 1) + " " + str(c + 1) + " " + str(d + 1) + " " + str(
                    e + 1) + " " + "0" + "\n"
                file.write(clauseseq)
                clauseseq = "-" + str(a + 1) + " " + str(b + 1) + " " + str(c + 1) + " " + str(d + 1) + " " + str(
                    e + 1) + " " + "0" + "\n"
                file.write(clauseseq)
                clauseseq = str(a + 1) + " " + str(b + 1) + " " + "-" + str(c + 1) + " " + "-" + str(
                    d + 1) + " " + "-" + str(e + 1) + " " + "0" + "\n"
                file.write(clauseseq)
                clauseseq = str(a + 1) + " " + "-" + str(b + 1) + " " + str(c + 1) + " " + "-" + str(
                    d + 1) + " " + "-" + str(e + 1) + " " + "0" + "\n"
                file.write(clauseseq)
                clauseseq = str(a + 1) + " " + "-" + str(b + 1) + " " + "-" + str(c + 1) + " " + str(
                    d + 1) + " " + "-" + str(e + 1) + " " + "0" + "\n"
                file.write(clauseseq)
                clauseseq = str(a + 1) + " " + "-" + str(b + 1) + " " + "-" + str(c + 1) + " " + "-" + str(
                    d + 1) + " " + str(e + 1) + " " + "0" + "\n"
                file.write(clauseseq)
                clauseseq = "-" + str(a + 1) + " " + str(b + 1) + " " + str(c + 1) + " " + "-" + str(
                    d + 1) + " " + "-" + str(e + 1) + " " + "0" + "\n"
                file.write(clauseseq)
                clauseseq = "-" + str(a + 1) + " " + str(b + 1) + " " + "-" + str(c + 1) + " " + str(
                    d + 1) + " " + "-" + str(e + 1) + " " + "0" + "\n"
                file.write(clauseseq)
                clauseseq = "-" + str(a + 1) + " " + str(b + 1) + " " + "-" + str(c + 1) + " " + "-" + str(
                    d + 1) + " " + str(e + 1) + " " + "0" + "\n"
                file.write(clauseseq)
                clauseseq = "-" + str(a + 1) + " " + "-" + str(b + 1) + " " + str(c + 1) + " " + str(
                    d + 1) + " " + "-" + str(e + 1) + " " + "0" + "\n"
                file.write(clauseseq)
                clauseseq = "-" + str(a + 1) + " " + "-" + str(b + 1) + " " + "-" + str(c + 1) + " " + str(
                    d + 1) + " " + str(e + 1) + " " + "0" + "\n"
                file.write(clauseseq)
                clauseseq = "-" + str(a + 1) + " " + "-" + str(b + 1) + " " + str(c + 1) + " " + "-" + str(
                    d + 1) + " " + str(e + 1) + " " + "0" + "\n"
                file.write(clauseseq)
                clauseseq = "-" + str(a + 1) + " " + "-" + str(b + 1) + " " + "-" + str(c + 1) + " " + "-" + str(
                    d + 1) + " " + "-" + str(e + 1) + " " + "0" + "\n"
                file.write(clauseseq)
        # 轮函数外的异或

        for k in range(2):
            for i in range(32):
                if k == 0:
                    a = xin[r][i]
                    b = xoutc[i]
                    c = xouta[i]
                else:
                    a = xin[r][i + 96]
                    b = xoutb[i]
                    c = xoutd[i]
                clauseseq = str(a + 1) + " " + str(b + 1) + " " + "-" + str(c + 1) + " " + "0" + "\n"
                file.write(clauseseq)
                clauseseq = str(a + 1) + " " + "-" + str(b + 1) + " " + str(c + 1) + " " + "0" + "\n"
                file.write(clauseseq)
                clauseseq = "-" + str(a + 1) + " " + str(b + 1) + " " + str(c + 1) + " " + "0" + "\n"
                file.write(clauseseq)
                clauseseq = "-" + str(a + 1) + " " + "-" + str(b + 1) + " " + "-" + str(c + 1) + " " + "0" + "\n"
                file.write(clauseseq)
    # Add constraints for the original sequential encoding
    Main_Vars = list([])
    for r in range(Round):
        for i in range(16):
            Main_Vars += [w[Round - 1 - r][i]]
    GenSequentialEncoding(Main_Vars, auxiliary_var_u, Main_Var_Num, CardinalityCons, file)
    # Add constraints for Matsui's strategy
    for matsui_count in range(0, MatsuiCount):
        StartingRound = MatsuiRoundIndex[matsui_count][0]
        EndingRound = MatsuiRoundIndex[matsui_count][1]
        LeftNode = 16 * StartingRound
        RightNode = 16 * EndingRound - 1
        PartialCardinalityCons = ActiveSbox - DiffActiveSbox[StartingRound] - DiffActiveSbox[Round - EndingRound]
        GenMatsuiConstraint(Main_Vars, auxiliary_var_u, Main_Var_Num, CardinalityCons, LeftNode, RightNode,
                            PartialCardinalityCons, file)
    file.close()

    # Call solver cadical
    order = "cadical " + "Problem-Round" + str(Round) + "-Active" + str(
        ActiveSbox) + ".cnf > Round" + str(Round) + "-Active" + str(ActiveSbox) + "-solution.out"
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
CountSbox = InitialLowerBound
TotalTimeStart = time.time()
for totalround in range(SearchRoundStart, SearchRoundEnd):
    flag = False
    time_start = time.time()
    MatsuiRoundIndex = []
    MatsuiCount = 0
    # Generate Matsui condition under choice 1
    if GroupConstraintChoice == 1:
        for group in range(0, GroupNumForChoice1):
            for round in range(1, totalround - group + 1):
                MatsuiRoundIndex.append([])
                MatsuiRoundIndex[MatsuiCount].append(group)
                MatsuiRoundIndex[MatsuiCount].append(group + round)
                MatsuiCount += 1
    # Printing Matsui conditions
    file = open("MatsuiCondition.out", "a")
    resultseq = "Round: " + str(totalround) + "; Partial Constraint Num: " + str(MatsuiCount) + "\n"
    file.write(resultseq)
    file.write(str(MatsuiRoundIndex) + "\n")
    file.close()
    # if totalround > 5:
    #     m = DiffActiveSbox[5]-DiffActiveSbox[3]  # branch
    #     roundm = DiffActiveSbox[totalround-2] + m
    #     if roundm>CountSbox:
    #         CountSbox = roundm
    # Flag2 = True
    # while Flag2:
    #     if MatsuiCount==0:
    #         Flag2=False
    #     for i in range(MatsuiCount):
    #         m1 = MatsuiRoundIndex[i][1]
    #         roundnumber = totalround - m1
    #         round_active_sbox = CountSbox - DiffActiveSbox[m1]
    #         if round_active_sbox < DiffActiveSbox[roundnumber]:
    #             CountSbox+=1
    #             break
    #         if i==MatsuiCount-1:
    #             Flag2=False

    while (flag == False):
        flag = Decision(totalround, CountSbox, MatsuiRoundIndex, MatsuiCount, flag)
        CountSbox += 1
    DiffActiveSbox[totalround] = CountSbox-1
    time_end = time.time()
    file = open("RunTimeSummarise.out", "a")
    resultseq = "Round: " + str(totalround) + "; Active S-box: " + str(
        DiffActiveSbox[totalround]) + "; Runtime: " + str(time_end - time_start) + "\n"
    print(time_end - TotalTimeStart)
    file.write(resultseq)
    file.close()
print(str(DiffActiveSbox))
TotalTimeEnd = time.time()
print("Total Runtime: " + str(TotalTimeEnd - TotalTimeStart))
file = open("RunTimeSummarise.out", "a")
resultseq = "Total Runtime: " + str(TotalTimeEnd - TotalTimeStart)
file.write(resultseq)
