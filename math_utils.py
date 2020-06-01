import numpy as np
import math

from itertools import combinations

def multiply_polynomials(greater_vec, less_vec):
    '''Multiplies two polynomials'''
    size = len(greater_vec)*len(less_vec)
    result = np.zeros(size, dtype=int)

    for i in range(len(greater_vec)):
        if greater_vec[i]:
            tmp_pol = np.concatenate((less_vec, np.zeros(len(greater_vec)-i-1, dtype=int)))
            result = xor(result, tmp_pol, size)

    return result

def multiply_bitwize(vec1, vec2, size: int):
    '''Multiplies two polynomials bitwize'''
    result = np.zeros(size, dtype=int)

    for i in range(size):
        result[i] = vec1[i] ^ vec2[i] 

    return result

def divide_polynomials(greater_vec, less_vec):
    '''Divides two polynomials, returns the quotient and the remainder'''
    res = 0
    len1 = len(greater_vec)
    quotient = np.zeros(len(greater_vec), dtype=int)
    remainder = greater_vec
    diff = greatest_bit(remainder) - greatest_bit(less_vec)

    while diff >= 0:
        tmp_pol = np.concatenate((less_vec, np.zeros(diff, dtype=int)))
        remainder = xor(remainder, tmp_pol, len(remainder))
        tmp_pol = np.zeros(diff + 1, dtype=int)
        tmp_pol[0] = 1
        quotient = xor(quotient, tmp_pol, len(remainder))
        diff = greatest_bit(remainder) - greatest_bit(less_vec)

    return quotient, remainder


def greatest_bit(vec):
    '''Calculates position of the greatest bit in vector'''
    for i in range(len(vec)):
        if vec[i]:
            return len(vec) - i
    return 0

def xor(vec1, vec2, k):
    '''Finds xor of two vectors'''
    result = np.zeros(k, dtype=int)

    vec1 = np.concatenate((np.zeros(k-len(vec1), dtype=int), vec1))
    vec2 = np.concatenate((np.zeros(k-len(vec2), dtype=int), vec2))

    for i in range(k):
        result[i] = vec1[i] ^ vec2[i]

    return result

def is_in(element, D):
    '''Determines if element is in given dictionary'''
    for vec in D:
        if element in vec:
            return True

    return False

def bin_vector_to_num(vec):
    '''Сonverts binary vector to number'''
    res = 0
    for i in range(len(vec)):
        res ^= (vec[i] << (len(vec) - i - 1))
    return res

def num_to_bin_vector(num):
    '''Converts number to binary vector'''
    if num == 0:
        return np.array([0])
    res = np.zeros(math.ceil(math.log2(num)) + 1, dtype = int)
    for i in range(math.ceil(math.log2(num)) + 1):
        if 1 & (num >> i):
            res[len(res) - i - 1] = 1
    return res[-greatest_bit(res):]

def get_powered_polynomial(polynomial, power):
    '''Finds the value of polynomial of the given argument power'''
    result = np.zeros(len(polynomial)*power, dtype=int)
    reversed_polynomial = polynomial[::-1]
    for i in range(len(reversed_polynomial)):
        if reversed_polynomial[i]:
            result[i*power] = 1
    result = result[::-1]
    if greatest_bit(result) == 0:
        return [0]
    else:
        return result[-greatest_bit(result):]