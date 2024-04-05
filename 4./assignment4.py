import sys
import time
import numpy as np


def get_input_data(filename):

    input_file = open(filename,'r')
    comment_list = list()
    gene_list = list()
    current_sequence = ''
    
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

    if current_sequence:
      gene_list.append(current_sequence)
    else :
       # 아무 데이터도 없는 경우
       gene_list = 'no DNA sequence'


    input_file.close()    

    # input data에 등장하는 모든 comment와 sequence list, return
    return gene_list





def edit_distance_measure(sequences):
    seq1 = sequences[0]
    seq2 = sequences[1] 

    m,n = len(seq1), len(seq2)

    if 2 != len(sequences):
        return 0
    else :
        dp = np.zeros((m + 1, n + 1), dtype=int)
    
    # 첫 번째 행과 열 초기화
    for i in range(m + 1):
        dp[i][0] = i
    for j in range(n + 1):
        dp[0][j] = j
    
    # use dynamic programming
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            cost = 0 if seq1[i - 1] == seq2[j - 1] else 1
            dp[i][j] = min(dp[i - 1][j] + 1,       # 삽입
                           dp[i][j - 1] + 1,       # 삭제
                           dp[i - 1][j - 1] + cost)  # 대체
    
    # 최소 편집 거리 반환
    return dp[m][n]


def output_data(filename, data):
    
    return 0



def main():
    input_filename = 'input.txt'
    output_filename = 'output.txt' 
    #input_filename = sys.argv[1]    
    #output_filename = sys.argv[2]

    dna_sequence = get_input_data(input_filename)

    #start_time= time.time()
    dp = edit_distance_measure(dna_sequence)
    #end_time = time.time()

    #elapsed_time = end_time - start_time

    print(dp)
    #output_data(output_filename,)
    
    #print(f"Elapsed Time: {elapsed_time:f} microseconds")


    return 0

if __name__ == '__main__':
    main()