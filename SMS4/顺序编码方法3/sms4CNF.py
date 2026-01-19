# 打开文件（替换 'your_file.txt' 为实际文件路径）
file_path = 'sms4.txt'
a= []
# 使用 with 语句打开文件，确保在处理完文件后自动关闭文件
with open(file_path, 'r', encoding='utf-8') as file:
    # 逐行读取文件内容
    for line in file:
        # 去除行尾的换行符（\n）
        tmp = []
        line = line.rstrip()
        for i in range(17):
            if line[i] == '-':
                tmp.append(9)
            elif line[i] == '1':
                tmp.append(1)
            else:
                tmp.append(0)
        # 处理每一行的内容
        print(tmp)
        a.append(tmp)
file2 = open('CNFConstraintForSbox.py', 'w')
file2.write('SymbolicCNFConstraintForSbox=')
file2.write(str(a))
print(a)


