# 주어진 4개의 문자열

# MSA_list -> 모두 global alignment로 완성된 sequence
def searching_same_char(MSA_list):

    strings = ["AsdD", "A-CD", "A-fD", "A---"]
    
    length = len(strings[0])
    star_list=list()
    
    for i in range(length):
        
        characters = [s[i] for s in strings]
        
        if all(char == characters[0] for char in characters):
            star_list.append("*")
        else:
            star_list.append(" ")


    for i in range(len(strings)):
        print(strings[i])

    result = ''.join(star_list)
    print(result)

    return result

