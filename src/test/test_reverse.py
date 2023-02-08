def reverse(number, length):
    reversed_num = 0;
    length >>= 1;

    while (length):
        reversed_num <<= 1;
        reversed_num |= number & 1;
        number >>= 1;
        length >>= 1;
    

    return reversed_num;

arr = [1, 1, 1, 2, 2, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0]
print(len(arr))
for i in range(16):
    print("i_reverse={}".format((reverse(i, 16))))
    print(arr[reverse(i, len(arr))])