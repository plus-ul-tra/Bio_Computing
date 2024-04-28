import numpy as np
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

    if first_char =='>':
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
    else:
        # FASTA format이 아닌 경우(#########했음#########)
        gene_list.append('nofasta')
        return gene_list

    if current_sequence:
      gene_list.append(current_sequence)
    else :
       # 아무 데이터도 없는 경우
       gene_list = 'no DNA sequence'


    input_file.close()    

    # input data에 등장하는 모든 comment와 sequence list, return
    return gene_list


# seqence, matching score 전달
def pairwise_sequence(seq1, seq2, match_score=2, mismatch_score=-2, gap_penalty=-1):
    n = len(seq1)
    m = len(seq2)

    # matrix 초기화
    score_matrix = np.zeros((n + 1, m + 1), dtype=int)

    # back trace 시 작성할 matrix 생성
    traceback_matrix = np.zeros((n + 1, m + 1), dtype=int)

    # 초기 gap penalty 설정
    for i in range(1, n + 1):
        score_matrix[i][0] = i * gap_penalty
        traceback_matrix[i][0] = 1
    for j in range(1, m + 1):
        score_matrix[0][j] = j * gap_penalty
        traceback_matrix[0][j] = 2

    # 각 자리간의 경우로 matrix 채우기
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match = score_matrix[i - 1][j - 1] + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_score)
            delete = score_matrix[i - 1][j] + gap_penalty
            insert = score_matrix[i][j - 1] + gap_penalty
            score_matrix[i][j] = max(match, delete, insert)
            if score_matrix[i][j] == match:
                # Diagonal
                traceback_matrix[i][j] = 0  

            elif score_matrix[i][j] == delete:
                # Up
                traceback_matrix[i][j] = 1  
            else:
                # Left
                traceback_matrix[i][j] = 2  

    # 최적경로 back trace
    aligned_seq1 = ''
    aligned_seq2 = ''
    i, j = n, m
    while i > 0 or j > 0:
        if i > 0 and j > 0 and traceback_matrix[i][j] == 0:
            aligned_seq1 = seq1[i - 1] + aligned_seq1
            aligned_seq2 = seq2[j - 1] + aligned_seq2
            i -= 1
            j -= 1
        elif i > 0 and traceback_matrix[i][j] == 1:
            aligned_seq1 = seq1[i - 1] + aligned_seq1
            aligned_seq2 = '-' + aligned_seq2
            i -= 1
        else:
            aligned_seq1 = '-' + aligned_seq1
            aligned_seq2 = seq2[j - 1] + aligned_seq2
            j -= 1

    return score_matrix[n][m], aligned_seq1, aligned_seq2

def output_data(filename, score, s1, s2):

    output_file =open(filename, 'w')

    if s1 =='empty':

        s='Empty file'
        output_file.write(s)
    elif s1 =='more':

        s='Need more sequences'
        output_file.write(s)

    elif s1 =='NoDNA':

        s='No DNA sequence'
        output_file.write(s)

    elif s1 =='nofasta':

        s='No FASTA format'
        output_file.write(s)

    else:
        s ='seqeunce1 : '+s1+'\n'+'sequence2 : '+ s2+'\n' + 'score :' + str(score)
        output_file.write(s)

    output_file.close()

    return 0

def main():
    #input_filename ='input.txt'
    #output_filename ='output.txt'
    input_filename = sys.argv[1]    
    output_filename = sys.argv[2]

    # command line argument로 input, output 파일 없을 시 예외 처리 (###########했음###########)
    if not os.path.exists(input_filename):
        print("No input file")
        return 0
        
    dna_sequence = get_input_data(input_filename)

    if len(dna_sequence) >= 2:
        # sequence 3 보다 많은 경우 처음 2개 사용

        # 대문자로 전부 치환 case insensitive 하도록
        seq1= dna_sequence[0].upper()
        seq2= dna_sequence[1].upper()

        # ATGC 가 아닌 다른 문자가 포함되는 경우 (###########했음###########)
        if any(base not in 'ATGC' for base in seq1) or any(base not in 'ATGC' for base in seq2):
            score =-100
            aligned_seq1 ='NoDNA'
            aligned_seq2 ='NoDNA'
        else:
            start_time= time.time()
            score, aligned_seq1, aligned_seq2 = pairwise_sequence(seq1, seq2)
            end_time = time.time()
            elapsed_time = (end_time - start_time) * 1000000
            # 정상 동작 시. 시간 출력(##########했음############)
            print(f"Elapsed Time: {elapsed_time:f} microseconds")
    #비정상 입력
    elif len(dna_sequence) ==1:
        # sequence 1개인 경우 score =-99, seq1, seq2에 'more' 대입 후 output 함수에서 처리
        # empty file의 경우 score = -100, seq1 , seq2 에 'empty' 대입 후 output 함수에서 처리

        # input이 빈 파일의 경우 (###########했음###########)
        if dna_sequence[0] =='file empty': 
            score = -100
            aligned_seq1 = 'empty'
            aligned_seq2 = 'empty'
        elif dna_sequence[0] =='nofasta':
            score=-102
            aligned_seq1 ='nofasta'
            aligned_seq2 ='nofasta'
        else:
            # dna sequence가 1개만 있을겨우 예외(###########했음###########)
            score = -99
            aligned_seq1 = 'more'
            aligned_seq2 = 'more'
    
    output_data(output_filename, score, aligned_seq1,aligned_seq2)


    return 0

if __name__ == '__main__':
    main()



