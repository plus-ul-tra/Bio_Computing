
import time
import sys
import os

def get_input_data(filename):
    input_file = open(filename,'r')
    comment_list = list()
    gene_list = list()
    current_sequence = ''
    
    first_char = input_file.read(1)
    input_file.seek(0)

    if first_char == '>':
    # \n를 고려하여 다음 comment, 즉 '>' 발견 전까지 모두 dna sequence로 설정.
        for line in input_file:
            line = line.strip()

            if line.startswith('>'):
                comment_list.append(line)

                if current_sequence:
                    gene_list.append(current_sequence)

                current_sequence = ''

            else:
                    current_sequence += line
    elif first_char =="":
        gene_list.append('nodata')
        code = -2
        return code, gene_list
    
    else:
        # FASTA format이 아닌 경우(#########했음#########)
        gene_list.append('nofasta')
        code = -1
        input_file.close()

        return code, gene_list

    if current_sequence:
      gene_list.append(current_sequence)
    else :
       # 아무 데이터도 없는 경우
       gene_list.append('nodata')
       code = -1
       return code, gene_list
       


    input_file.close()    

    # 예외코드 -3, -4는 해당 함수에서 선별
    exception_code, final_data = exception_check(gene_list)

    # input data에 등장하는 모든 comment와 sequence list, return
    return exception_code, final_data




def exception_check(data_list):
    # 여기까지 넘어 왔다. fasta format을 만족하는 데이터 들이 있긴하다
    # 동작 : list 개수, upper 로 변환,
    new_data_list = list()
    code = 0
    if len(data_list) == 1:
        # data가 하나 뿐이다 최소 2개는 있어야함
        code = -3
        return code, data_list
    
    for data in data_list:
        new_data = data.upper()
        if any(base not in 'ATGC' for base in new_data):
            code = -4

            return code, data_list
        else:
            new_data_list.append(new_data)
    
        
    # 아무 예외 없는 정상 데이터의 경우 code = 0
    return code, new_data_list

# 정상 알고리즘 출력
def result_output(filename, data,length):
    output_file = open(filename, 'w')
    if data:
        s = 'length '+ str(length) +' pattern : ' + str(data)

    else:
        s = 'No pattern found'

    output_file.write(s)
    output_file.close()
    return 0



def exception_output(filename, code):
    output_file = open(filename, 'w')
    # -1 : no fasta format
    # -2 : no data in file
    # -3 : single sequence
    # -4 : no dna sequence
    if code == -1:
        s = 'No correct format'
    elif code == -2:
        s = 'No DNA sequence'
    elif code == -3:    
        s = 'Need more sequences'
    elif code == -4:
        s = 'No DNA sequence'

    output_file.write(s)
    output_file.close()    

    return 0


def approximate_pattern_finding(sequences,pattern_length):

    subsequences_list = [extract_subsequences(seq, pattern_length) for seq in sequences]

    min_distance = float('inf')
    best_pattern = set()

    for i in range(len(subsequences_list)):
        for j in range(i + 1, len(subsequences_list)):
            for subseq1 in subsequences_list[i]:
                for subseq2 in subsequences_list[j]:
                    distance = levenshtein_distance(subseq1, subseq2)
                    if distance < min_distance:
                        min_distance = distance
                        best_pattern.add(subseq1)
                        best_pattern.add(subseq2)
    print(best_pattern)
    return best_pattern

def levenshtein_distance(s1, s2):
    len_s1 = len(s1)
    len_s2 = len(s2)

    dp = [[0] * (len_s2 + 1) for _ in range(len_s1 + 1)]

    for i in range(len_s1 + 1):
        dp[i][0] = i
    for j in range(len_s2 + 1):
        dp[0][j] = j

    for i in range(1, len_s1 + 1):
        for j in range(1, len_s2 + 1):
            if s1[i - 1] == s2[j - 1]:
                dp[i][j] = dp[i - 1][j - 1]
            else:
                dp[i][j] = min(dp[i - 1][j] + 1,
                               dp[i][j - 1] + 1,
                               dp[i - 1][j - 1] + 1)

    return dp[len_s1][len_s2]

def extract_subsequences(sequence, length):
    subsequences = []
    seen = set()  

    for i in range(len(sequence) - length + 1):
        subsequence = sequence[i:i + length]
        # 만약 서브시퀀스가 이미 존재하는지 확인하고, 존재하지 않으면 추가
        if subsequence not in seen:
            subsequences.append(subsequence)
            seen.add(subsequence)

    return subsequences


def main():
    #input_filename ='input.txt'
    #output_filename ='output.txt'
    #pattern_length = 8
    input_filename = sys.argv[1]    
    output_filename = sys.argv[2]
    pattern_length = int(sys.argv[3])
    if not os.path.exists(input_filename):
        print("No input file")
        return 0
    
    ex_code, dna_sequence = get_input_data(input_filename)

    if ex_code != 0:
        exception_output(output_filename,ex_code)
        return 0
    start_time= time.time()
    result = approximate_pattern_finding(dna_sequence, pattern_length)
    end_time = time.time()
    elapsed_time = (end_time - start_time) * 1000000
    print(f"Elapsed Time: {elapsed_time:f} microseconds")
    result_output(output_filename, result, pattern_length)
    
    return 0




if __name__ == '__main__':
    main()