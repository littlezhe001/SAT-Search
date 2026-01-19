import copy
import os
import time

from CNFConstraintForSbox import SymbolicCNFConstraintForSbox
FullRound = 30

SearchRoundStart = 1
SearchRoundEnd = 30
InitialLowerBound =1
GroupConstraintChoice = 1

# Parameters for choice 1
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
        # count +=64*24+32*36
        count +=8286*4 + 32*28

    return count


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


def genclause(X, i, l1):
    ibin = bin(i)[2:]
    while len(ibin) != l1:
        ibin = '0' + ibin
    clause = ''
    for i in range(l1):
        if ibin[i] == '0':
            clause = clause + str(X[i] + 1) + ' '
        else:
            clause = clause + '-' + str(X[i] + 1) + ' '
    w = 1
    for i in range(l1):
        w = w ^ int(ibin[i], 2)
    if w:
        clause = clause + '-' + str(X[-1] + 1) + ' '
    else:
        clause = clause + str(X[-1] + 1) + ' '
    clause = clause + '0' + '\n'
    return clause




def Decision(Round, ActiveSbox, MatsuiRoundIndex, MatsuiCount, flag):
    TotalSbox = 4 * Round
    count_var_num = 0
    time_start = time.time()
    xin = []
    xout = []
    sin = []
    sout = []
    w = []
    xor1 = []
    xor2 = []
    for i in range(Round):
        xin.append([])
        xout.append([])
        sin.append([])
        sout.append([])
        xor1.append([])
        xor2.append([])
        w.append([])
        for j in range(128):
            xin[i].append(0)
            xout[i].append(0)
        for j in range(32):
            sin[i].append(0)
            sout[i].append(0)
            xor1[i].append(0)
            xor2[i].append(0)
        for j in range(4):
            w[i].append(0)

    # Allocate variables
    for j in range(128):
        xin[0][j] = count_var_num
        count_var_num += 1
    for i in range(1,Round):
        for j in range(96):
            xin[i][j] = xin[i-1][j+32]
        for j in range(96,128):
            xin[i][j] = count_var_num
            count_var_num += 1
    for i in range(Round):
        for j in range(32):
            sin[i][j] = count_var_num
            count_var_num += 1
        for j in range(32):
            sout[i][j] = count_var_num
            count_var_num += 1
        for j in range(32):
            xor1[i][j] = count_var_num
            count_var_num +=1
        for j in range(32):
            xor2[i][j] = count_var_num
            count_var_num +=1
        for j in range(4):
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
        if i < ActiveSbox:
            for j in range(i + 1):
                auxiliary_var_u[i].append(count_var_num)
                count_var_num += 1
        else:
            for j in range(ActiveSbox):
                auxiliary_var_u[i].append(count_var_num)
                count_var_num += 1

    # Count the number of clauses in the round function
    count_clause_num = 0
    count_clause_num = CountClausesInRoundFunction(Round, ActiveSbox, count_clause_num)
    # Count the number of clauses in the original sequential encoding
    Main_Var_Num = 4 * Round
    CardinalityCons = ActiveSbox
    count_clause_num = CountClausesInSequentialEncoding(Main_Var_Num, CardinalityCons, count_clause_num)
    # Count the number of clauses for Matsui's strategy
    for matsui_count in range(0, MatsuiCount):
        StartingRound = MatsuiRoundIndex[matsui_count][0]
        EndingRound = MatsuiRoundIndex[matsui_count][1]
        LeftNode = 4 * StartingRound
        RightNode = 4 * EndingRound - 1
        PartialCardinalityCons = ActiveSbox - DiffActiveSbox[StartingRound] - DiffActiveSbox[Round - EndingRound]
        count_clause_num = CountClausesForMatsuiStrategy(Main_Var_Num, CardinalityCons, LeftNode, RightNode,
                                                         PartialCardinalityCons, count_clause_num)

    # Open file
    file = open("Problem-Round" + str(Round) + "-Active" + str(ActiveSbox) + ".cnf", "w")
    file.write("p cnf " + str(count_var_num) + " " + str(count_clause_num) + "\n")
    print("所需变量个数为：" + str(count_var_num) + " " + "所需要的句子数量为：" + str(count_clause_num))
    # Add constraints to claim nonzero input difference
    clauseseq = ""
    for i in range(128):
        clauseseq += str(xin[0][i] + 1) + " "
    clauseseq += "0" + "\n"
    file.write(clauseseq)

    # Add constraints for the round function

    for r in range(Round):
        # 进入S盒前的异或
        x0 = xin[r][:32]
        x1 = xin[r][32:64]
        x2 = xin[r][64:96]
        x3 = xin[r][96:128]
        inputxor1 = sin[r]
        for i in range(32):
            X = [x1[i],x2[i],x3[i],inputxor1[i]]
            for k in range(8):
                clauseseq = genclause(X,k,3)
                file.write(clauseseq)
        # S盒约束
        for i in range(4):
            X = list([])
            for k in range(8):
                X += [sin[r][8 * i + k]]
            for k in range(8):
                X += [sout[r][8 * i + k]]
            X += [w[r][i]]
            for j in range(8286):
                clauseseq = ""
                for k in range(17):
                    if (SymbolicCNFConstraintForSbox[j][k] == 1):
                        clauseseq += "-" + str(X[k] + 1) + " "
                    if (SymbolicCNFConstraintForSbox[j][k] == 0):
                        clauseseq += str(X[k] + 1) + " "
                clauseseq += "0" + "\n"
                file.write(clauseseq)

        # L变换约束
        y0 = sout[r]
        y1 = sout[r][2:]+sout[r][:2]
        y2 = sout[r][10:]+sout[r][:10]
        y3 = sout[r][18:]+sout[r][:18]
        y4 = sout[r][24:]+sout[r][:24]
        inputxor2 = xor1[r]
        inputxor3 = xor2[r]

        for i in range(32):
            X = [y0[i],y1[i],y2[i],inputxor2[i]]
            for k in range(8):
                clauseseq = genclause(X,k,3)
                file.write(clauseseq)
        for i in range(32):
            X = [inputxor2[i],y3[i],y4[i],inputxor3[i]]
            for k in range(8):
                clauseseq = genclause(X,k,3)
                file.write(clauseseq)

        inputxor4 = xout[r][96:]
        for i in range(32):
            X = [x0[i],inputxor3[i],inputxor4[i]]
            for k in range(4):
                clauseseq = genclause(X,k,2)
                file.write(clauseseq)



    # Add constraints for the original sequential encoding
    Main_Vars = list([])
    for r in range(Round):
        for i in range(4):
            Main_Vars += [w[Round - 1 - r][i]]
    GenSequentialEncoding(Main_Vars, auxiliary_var_u, Main_Var_Num, CardinalityCons, file)
    # Add constraints for Matsui's strategy
    for matsui_count in range(0, MatsuiCount):
        StartingRound = MatsuiRoundIndex[matsui_count][0]
        EndingRound = MatsuiRoundIndex[matsui_count][1]
        LeftNode = 4 * StartingRound
        RightNode = 4 * EndingRound - 1
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
            for round in range(0, totalround - group - 1):
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
    # 外面也限制
    # Flag2 = True
    # while Flag2:
    #     if MatsuiCount == 0:
    #         Flag2 = False
    #     for i in range(MatsuiCount):
    #         m1 = MatsuiRoundIndex[i][1]
    #         roundnumber = totalround - m1
    #         round_active_sbox = CountSbox - DiffActiveSbox[m1]
    #         if round_active_sbox < DiffActiveSbox[roundnumber]:
    #             CountSbox += 1
    #             break
    #         if i == MatsuiCount - 1:
    #             Flag2 = False

    while (flag == False):
        flag = Decision(totalround, CountSbox, MatsuiRoundIndex, MatsuiCount, flag)
        CountSbox += 1
    CountSbox -=1
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
