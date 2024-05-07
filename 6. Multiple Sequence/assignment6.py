import numpy as np
import time
import sys
import os

def get_blosum62_matrix(filename):
    blosum62 = open(filename, 'r')
    amino_list = list()
    score_matrix = list()
    lines = blosum62.readlines()[1:]
    for line in lines:
        row=list(line[1:].split())
        score_matrix.append(row)

    return score_matrix


# cmd + shift + L
def get_input_data(filename):
    input_file = open(filename,'r')
    comment_list = list()
    protein_list = list()
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
                    protein_list.append(current_sequence)

                current_sequence = ''

            else:
                    current_sequence += line
    else:
        # FASTA format이 아닌 경우
        protein_list.append('nofasta')
        return protein_list

    if current_sequence:
      protein_list.append(current_sequence)
    else :
       # 아무 데이터도 없는 경우
       protein_list = 'no DNA sequence'


    input_file.close()    

    # input data에 등장하는 모든 comment와 sequence list, return
    return protein_list


# pairwise alignment 쌍에 대한 blosum62를 이용하여 scoring
def get_blosum_score(amino1, amino2,matrix):
    amino_acid = "CSTPAGNDEQHRKMILVFYW"
    index1 = amino_acid.find(amino1)
    index2 = amino_acid.find(amino2)
    score = int(matrix[index1][index2])
    
    return score


def pairwise_sequence(seq1, seq2, matrix, gap_penalty = -5):
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
            # match 된 경우에서 뒤에 더해지는 score를 blosum62의 table을 매김.
            # 이부분만 수정 -> seq1[i-1], seq2[j-1]을 받아서 blosum62의 점수 반환하는 함수사용
            value = get_blosum_score(seq1[i - 1],seq2[j - 1],matrix)
            # (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_score)
            match = score_matrix[i - 1][j - 1] + value
            # delete, insert는 모두 gap (옆, 밑 이동을 의미)
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
    print(score_matrix)
    print(traceback_matrix)
    return score_matrix[n][m], aligned_seq1, aligned_seq2

def select_center_sequence(sequence_list, matrix):

    max_strength = 0
    for s1 in sequence_list:
        strength = 0
        for s2 in sequence_list:
            if s1 != s2:
                #본인을 제외한 나머지 sequence와 pairwise 구함.
                score, align1,align2 = pairwise_sequence(s1,s2, matrix) 
                strength += score
        
        if strength > max_strength :
            max_strength = strength
            center_sequence = s1


    return center_sequence

def main():
    input_filename ='input.txt'
    output_filename ='output.txt'
    blosum62_file = 'blosum62.txt'
    #input_filename = sys.argv[1]    
    #output_filename = sys.argv[2]
    #blosum62 생성
    matrix = get_blosum62_matrix(blosum62_file)
    ex_list = list()
    
    # 점수가 높을수록 유사도가 깊음
    s1 = "CSTPGVKWCCKWCTT"
    s2 = "CTTPGVVKCCKRT"
    s3 = "TTCGVVWkTT"

    ex_list.append(s1)
    ex_list.append(s2)
    ex_list.append(s3)
    #score, r1,r2 = pairwise_sequence(s1,s3,matrix)
    center_sequence = select_center_sequence(ex_list, matrix)

    return 0

if __name__ == '__main__':
    main()