import numpy as np
import math_utils
from itertools import combinations


class GaloisField:
    '''Represents Galois field GF(2^p) and related operations
       
       Parameters:
        p - determines the power of the field 
     '''

    def __init__(self, p: float):
        self.field_power = 2 ** p -1
        self.power = p
        self.cyclotomic_classes = self.get_cyclotomic_classes()
        self.primitive_polynomial = self.get_primitive_polynomial()
        self.logarithm_table = self.get_logarithm_table()
        self.logarithm_table[-1] = [0]
        self.reversed_log_table = self.reverse_log_table()
        self.reversed_log_table[0] = -1


    def get_primitive_polynomial(self):
        primitive_polynomial = {
            2 : 0b111,
            3 : 0b1011,
            4 : 0b10011,
            5 : 0b100101,
            6 : 0b1000011,
            7 : 0b10001001,
            8 : 0b100011101,
            9 : 0b1000010001,
            10 : 0b10000001001,
            11 : 0b100000000101,
            12 : 0b1000001010011,
            13 : 0b10000000011011,
            14 : 0b100010001000011,
            15 : 0b1000000000000011,
            16 : 0b10001000000001011,
            17 : 0b100000000000001001,
            18 : 0b1000000000010000001,
            19 : 0b10000000000000100111,
            20 : 0b100000000000000001001
        }
        return np.array(list(np.binary_repr(primitive_polynomial[self.power], self.power+1)), dtype=int)


    def get_cyclotomic_classes(self):
        '''Calculates the cyclotomic classes in the given field'''
        classes = {}
        num = 0

        for i in range(self.field_power):
            if math_utils.is_in(i, classes.values()):
                continue

            classes[num] = [i]
            tmp = i * 2 % (self.field_power)
            while tmp not in classes[num]:
                classes[num].append(tmp)
                tmp = tmp * 2 % (self.field_power)
            num += 1

        return classes


    def get_logarithm_table(self):
        '''Calculates powers of the primitive element in the given field'''
        size = math_utils.greatest_bit(self.primitive_polynomial) - 1
        logarithm_table = {}

        for i in range(size):
            logarithm_table[i] = np.zeros(i+1, dtype=int)
            logarithm_table[i][0] = 1

        logarithm_table[size] = self.primitive_polynomial[1:]

        for i in range(size + 1, 2 ** size - 1):
            tmp_pol = np.concatenate((logarithm_table[i-1], np.zeros(1, dtype=int)))
            if tmp_pol[len(tmp_pol)-size-1]:
                tmp_pol = math_utils.xor(tmp_pol, self.primitive_polynomial[1:], size+1)
            logarithm_table[i] = tmp_pol[1:]
        return logarithm_table


    def get_minimal_polynomial(self, cycl_class):
        '''Calculates the minimal polynomial for the given cyclotomic class'''
        polynomial = 1 << len(cycl_class)
        polynomial = np.concatenate(([1], np.zeros(len(cycl_class), dtype = int)))
        for i in range(len(cycl_class)):
            coefficient = np.zeros(self.power, dtype = int)
            for combination in combinations(cycl_class, i + 1):
                coefficient = math_utils.xor(coefficient, 
                        self.logarithm_table[sum(combination) % self.field_power], self.power)
            if (sum(coefficient) != 0):
                coefficient = coefficient[-math_utils.greatest_bit(coefficient):]
                coefficient = np.concatenate((coefficient, np.zeros(len(cycl_class) - i - 1, dtype = int)))
                polynomial = math_utils.xor(coefficient, polynomial, len(cycl_class) + 1)

        return polynomial


    def find_roots_of_polynomial(self, polynomial):
        '''Calculates the roots of the given polynomial'''
        roots = []
        for candidate in range(self.field_power):
            result = self.logarithm_table[polynomial[0]]
            for polynomial_power in range(1, len(polynomial)):
                if polynomial[polynomial_power] >= 0:
                    temp = self.logarithm_table[(polynomial[polynomial_power] + candidate * polynomial_power) % (self.field_power)]
                    result = math_utils.xor(temp, result, max(len(result), len(temp)))
            if sum(result) == 0:
                roots.append(candidate)

        return roots

    def reverse_log_table(self):
        '''Reverses the logarithm table so that the key of the table
           is the value of the power'''
        new_table = dict()
        for key, value in self.logarithm_table.items():
            new_key = math_utils.bin_vector_to_num(value)
            new_table[new_key] = key

        return new_table
        